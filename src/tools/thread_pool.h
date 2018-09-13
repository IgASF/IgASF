/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once
// https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::result_of<F(Args...)>::type>;
    ~ThreadPool();
    void waitAll();
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void()> > tasks;
    // synchronization
    std::mutex wait_mutex;
    std::mutex queue_mutex;
    std::condition_variable condition;
    std::condition_variable wait_var;
    bool stop;
    int running;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    :   stop(false), running(0)
{
    for(size_t i = 0;i<threads;++i)
        workers.emplace_back(
            [this]
            {
                while (!this->stop)
                {
                    std::function<void()> task;
//                    {
//                    std::unique_lock<std::mutex> lock(this->queue_mutex);
//                    if(tasks.empty()) {
//                        this->condition.wait_for(lock, std::chrono::duration<int, std::milli>(10));
//                        continue;
//                    }
//                    task = std::move(this->tasks.front());
//                    this->tasks.pop();
//                    }
                
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });
                        if(this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    {
                        std::unique_lock<std::mutex> lock(this->wait_mutex);
                        this->running+=1;
                    }
                    task();
                    {
                        std::unique_lock<std::mutex> lock(this->wait_mutex);
                        this->running-=1;
                    }
                    this->wait_var.notify_one();
                }
            }
        );
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::future<typename std::result_of<F(Args...)>::type>
{
    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared< std::packaged_task<return_type()> >(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        // don't allow enqueueing after stopping the pool
        if(stop)  throw std::runtime_error("enqueue on stopped ThreadPool");
        tasks.emplace([task](){ (*task)(); });
    }
    condition.notify_one();
    return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

void ThreadPool::waitAll()
{
    std::unique_lock<std::mutex> lock(this->wait_mutex);
    this->wait_var.wait(lock, [this]{ return this->running==0 && this->tasks.empty();; });
}

#endif
