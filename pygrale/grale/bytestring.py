import sys

# Set up functions B and S to convert a string to bytes or vice versa.
B = lambda s: bytes(s, 'UTF-8')
S = lambda b: b.decode(encoding='UTF-8')


