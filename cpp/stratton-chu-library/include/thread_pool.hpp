#pragma once

#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <vector>
#include <algorithm>
#include "function.hpp"

template <typename T>
class ConcurentQueue {
    std::mutex m_mut;
    std::queue<T> m_q;

public:
    using value_type = T;

    ConcurentQueue() = default;
    ConcurentQueue(const ConcurentQueue& other) : m_mut{}, m_q{other.m_q} {}

    void push(T el) {
        std::lock_guard lk{m_mut};
        m_q.push(el);
    }

    T pop() {
        std::lock_guard lk{m_mut};
        if (m_q.size() > 0) {
            auto ret = m_q.front();
            m_q.pop();
            return ret;
        }
        else {
            return T{};
        }
    }
};

template <typename T>
static void worker(ConcurentQueue<T>& task_queue){
    while (true) {
        auto task = task_queue.pop();
        if (task) task();
        else break;
    }
};


template <typename T>
class ThreadPoolImpl {
    T m_task_queue;

public:
    ThreadPoolImpl(const T& task_queue)
        : m_task_queue{task_queue}
    {}

    template <typename QueueItem>
    void run(void(*worker)(ConcurentQueue<QueueItem>&)) {
        const auto system_thread_count = std::thread::hardware_concurrency();
        std::vector<std::thread> pool{system_thread_count - 2u};
        std::for_each(pool.begin(), pool.end(),
            [&](auto& thread){ thread = std::thread{worker, std::ref(m_task_queue)}; });
        std::for_each(pool.begin(), pool.end(), [&](auto& thread){ thread.join(); });
    }
};
