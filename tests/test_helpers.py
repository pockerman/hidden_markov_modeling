"""
Unit tests for helpers module
"""
import unittest
import json

from helpers import partition_range


class TestHelpers(unittest.TestCase):


  def test_equal_partition(self):
    start = 0
    end = 10
    pieces = 2
    chuncks = partition_range(start=start, end=end, npieces=pieces)

    self.assertEqual(len(chuncks), 2)
    self.assertEqual(chuncks[0][0], 0)
    self.assertEqual(chuncks[0][1], 4)
    self.assertEqual(chuncks[1][0], 5)
    self.assertEqual(chuncks[1][1], 10)


  def test_unequal_partition(self):
    start = 0
    end = 15
    pieces = 4
    chuncks = partition_range(start=start, end=end, npieces=pieces)

    self.assertEqual(len(chuncks), 4)
    self.assertEqual(chuncks[0][0], 0)
    self.assertEqual(chuncks[0][1], 2)
    self.assertEqual(chuncks[1][0], 3)
    self.assertEqual(chuncks[1][1], 5)
    self.assertEqual(chuncks[2][0], 6)
    self.assertEqual(chuncks[2][1], 8)
    self.assertEqual(chuncks[3][0], 9)
    self.assertEqual(chuncks[3][1], 15)




if __name__ == '__main__':
    unittest.main()



