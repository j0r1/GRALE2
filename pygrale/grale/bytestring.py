import sys

# Set up functions B and S to convert a string to bytes or vice versa. Useful
# for both python 2 and 3 support
if sys.version_info[0] == 2:
    B = lambda s: s
    S = B
else:
    B = lambda s: bytes(s, 'UTF-8')
    S = lambda b: b.decode(encoding='UTF-8')


