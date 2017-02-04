#ifndef CONCURRENT_QUEUE_H
#define CONCURRENT_QUEUE_H
#include <queue>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>
#include <vector>
template<typename Data>
class concurrent_queue
{
private:
	std::queue<Data> the_queue;
	mutable boost::mutex the_mutex;
	boost::condition_variable the_condition_variable;
	std::vector<int> batch;
	std::vector<int> single;
public:

	void set_batch_number(int amount) {
		for(int i = 0; i < amount; i++) {
			batch.push_back(0);
		}
	}

	void set_single_number(int amount) {
		for(int i = 0; i < amount; i++) {
			single.push_back(0);
		}
	}

	std::vector<int>* get_single_set() {
		return &single;
	}

	std::vector<int>* get_batch_set() {
		return &batch;
	}

	void push(Data const& data)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		the_queue.push(data);
		lock.unlock();
		the_condition_variable.notify_one();
	}

	bool empty() const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		return the_queue.empty();
	}

	int size() const
	{
		boost::mutex::scoped_lock lock(the_mutex);
		return the_queue.size();
	}

	bool try_pop(Data& popped_value, int id)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return false;
		}
		single[id] += 1;
		popped_value=the_queue.front();
		the_queue.pop();
		return true;
	}

	int try_pop_batch(std::vector<int>* data, int amount, int id)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		if(the_queue.empty())
		{
			return 0;
		}
		int i = 0;
		for( ; i < amount && !the_queue.empty(); i++){
			int popped_value= the_queue.front();
			the_queue.pop();
			data->push_back(popped_value);
		}
		batch[id] += i;
		return i;
	}

	void wait_and_pop(Data& popped_value)
	{
		boost::mutex::scoped_lock lock(the_mutex);
		while(the_queue.empty())
		{
			the_condition_variable.wait(lock);
		}

		popped_value=the_queue.front();
		the_queue.pop();
	}

	void clear() {

		boost::mutex::scoped_lock lock(the_mutex);
		std::queue<Data> empty;
		std::swap( the_queue, empty );
		lock.unlock();
	}

};
#endif // CONCURRENT_QUEUE_H

