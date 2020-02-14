"""
python -m tests.exports.test_task_scheduler
"""
import unittest
import random

import ee

from gee_tools.exports.task_scheduler import TaskScheduler

# This may need to inherit from Task
class FakeTask(object):

    def __init__(self, time_until_finished, fail=False):
        self.started = False
        self.remaining = time_until_finished
        self.fail = fail

    def start(self):
        if self.started:
            raise RuntimeError('Task already started.')
        self.started = True

    def _status(self):
        if not self.started:
            return ee.batch.Task.State.UNSUBMITTED
        if self.remaining > 0:
            self.remaining -= 1
            return ee.batch.Task.State.RUNNING
        if self.fail:
            return ee.batch.Task.State.FAILED
        return ee.batch.Task.State.COMPLETED

    def status(self):
        return {'state': self._status()}

    @property
    def State(self):
        return ee.batch.Task.State


def add_task(scheduler, jid, dependencies):
    dependencies = list(dependencies)
    num_dep = 0 if len(dependencies) == 0 else random.randint(0, (len(dependencies) - 1))
    random.shuffle(dependencies)
    dependencies = dependencies[:num_dep]

    fail = bool(random.randint(0, 1))
    time_until_finished = random.randint(1, 10)

    task = FakeTask(time_until_finished, fail=fail)
    scheduler.add_task(task, jid, dependencies)


def test(n):
    random.seed(0)
    scheduler = TaskScheduler()
    for i in range(n):
        add_task(scheduler, i, list(range(i)))
    if bool(random.randint(0, 1)):
        scheduler.mark_completed(n // 2)
    scheduler.run(max_processes=8, sleep_time=0.0, verbose=0)


class TaskSchedulerUnitTest(unittest.TestCase):

    def test(self):
        for _ in range(10):
            test(10)
        for _ in range(50):
            test(100)

    def test_exactly_4(self):
        random.seed(0)
        scheduler = TaskScheduler()
        for i in range(4):
            task = FakeTask(10, fail=False)
            scheduler.add_task(task, i)
        scheduler.run(max_processes=4, sleep_time=0.0, verbose=0)


if __name__ == '__main__':
    unittest.main()
