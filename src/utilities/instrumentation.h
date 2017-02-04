#ifndef INSTRUMENTATION_H
#define INSTRUMENTATION_H

struct timer {
    timespec time_start;
    timespec time_end;

    void start() {
        clock_gettime(CLOCK_MONOTONIC, &time_start);
    }

    void stop() {
        clock_gettime(CLOCK_MONOTONIC, &time_end);
    }

    double elapsed() {
        return (time_end.tv_sec - time_start.tv_sec)
                + (time_end.tv_nsec	- time_start.tv_nsec) * 1e-9;
    }

};

#endif // INSTRUMENTATION_H

