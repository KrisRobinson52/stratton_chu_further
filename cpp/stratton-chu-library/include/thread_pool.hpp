#pragma once

#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <vector>
#include <algorithm>

template <typename T>
class ConcurentQueue {
    std::mutex m_mut;
    std::queue<T> m_q;

public:
    ConcurentQueue() = default;
    ConcurentQueue(const ConcurentQueue& other) : m_mut{}, m_q{other.m_q} {}

    void push(T el) {
        std::scoped_lock lk(m_mut);
        m_q.push(el);
    }

    T pop() {
        std::scoped_lock lk(m_mut);
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
class ThreadPoolImpl {
    T m_task_queue;
    std::vector<std::thread> m_threads;

public:
    ThreadPoolImpl(const T& task_queue)
        : m_task_queue{task_queue},
        m_threads{std::thread::hardware_concurrency() - 2}  // Устанавлиевает количество потоков как количество системных потоков - 2
        // m_threads{1}  // В один поток
    {}

    void run() {
        const auto task = [&](){
            while (true) {
                auto task = m_task_queue.pop();
                if (task) task();
                else break;
            }
        };
        std::for_each(std::begin(m_threads), std::end(m_threads),
            [&](auto& thread){ thread = std::thread{task}; });
        std::for_each(std::begin(m_threads), std::end(m_threads), [&](auto& thread){ thread.join(); });
    }
};


using TaskQueue = ConcurentQueue<std::function<void()>>;
using ThreadPool = ThreadPoolImpl<ConcurentQueue<std::function<void()>>>;
