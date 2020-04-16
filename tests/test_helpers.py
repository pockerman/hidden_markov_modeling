"""
Unit tests for helpers module
"""
import unittest
import json

from helpers import Observation
from helpers import observation_to_json
from helpers import Window
from helpers import windows_to_json
from helpers import windows_from_json
from helpers import listify_dicts_property


class TestHelpers(unittest.TestCase):


  def test_observarion_to_json(self):
    """
    Test if an observation is dumped correctly
    into a json formatted string

    Returns
    -------
    None.

    """

    observation = Observation(position=100,
                              read_depth=200,
                              base=['A','B'])

    json_str = observation_to_json(observation)
    data = json.loads(json_str)

    self.assertEqual(data["position"], 100)
    self.assertEqual(data["read_depth"], 200)
    self.assertEqual(data["base"], ['A','B'])


  def test_window_to_json(self):
    """
    Test if a window is properly dumped into a
    json formatted string

    Returns
    -------
    None.

    """

    window = Window(idx=10, capacity=5)
    window.add(observation=Observation(position=100,
                              read_depth=200,
                              base=['A','B']))

    window.add(observation=Observation(position=300,
                              read_depth=500,
                              base=['A','B']))

    json_str = window.to_json()
    data = json.loads(json_str)

    self.assertEqual(data["id"], window.get_id())
    self.assertEqual(data["capacity"], window.capacity())
    self.assertEqual(len(data["observations"]), 2)
    self.assertEqual(data["observations"][0],
                     '{"position": 100, "read_depth": 200, "base": ["A", "B"]}')
    self.assertEqual(data["observations"][1],
                     '{"position": 300, "read_depth": 500, "base": ["A", "B"]}')

  def test_windows_to_json(self):

    windows = []
    window = Window(idx=10, capacity=5)
    window.add(observation=Observation(position=100,
                              read_depth=200,
                              base=['A','B']))

    window.add(observation=Observation(position=300,
                              read_depth=500,
                              base=['A','B']))

    windows.append(window)

    window = Window(idx=11, capacity=5)
    window.add(observation=Observation(position=100,
                              read_depth=200,
                              base=['A','B']))

    window.add(observation=Observation(position=300,
                              read_depth=500,
                              base=['A','B']))

    windows.append(window)
    json_str = windows_to_json(windows=windows)
    data = json.loads(json_str)
    new_windows = windows_from_json(jsonmap=data)

    self.assertEqual(new_windows[0].get_id(), windows[0].get_id())
    self.assertEqual(new_windows[0].capacity(), windows[0].capacity())
    self.assertEqual(new_windows[0][0].position, windows[0][0].position)
    self.assertEqual(new_windows[0][0].read_depth, windows[0][0].read_depth)
    self.assertEqual(new_windows[0][0].base, windows[0][0].base)

    self.assertEqual(new_windows[1].get_id(), windows[1].get_id())
    self.assertEqual(new_windows[1].capacity(), windows[1].capacity())
    self.assertEqual(new_windows[1][0].position, windows[1][0].position)
    self.assertEqual(new_windows[1][0].read_depth, windows[1][0].read_depth)
    self.assertEqual(new_windows[1][0].base, windows[1][0].base)

  def test_listify_dicts_property(self):

    test_list = [{"prop1":100, "prop2":200},
                 {"prop1":500, "prop3":200}]


    values = listify_dicts_property(list_dict_vals=test_list,
                                    property_name="prop1")


    self.assertEqual(len(values), 2)
    self.assertEqual(values[0], 100)
    self.assertEqual(values[1], 500)


if __name__ == '__main__':
    unittest.main()



