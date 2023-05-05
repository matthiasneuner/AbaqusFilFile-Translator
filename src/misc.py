from collections import defaultdict


class RecursiveDefaultDict(dict):
    # return defaultdict(RecursiveDefaultDict)

    def __init__(self, maxLevels):
        self.level = maxLevels
        return super().__init__()

    def __getitem__(self, key):
        if key in self:
            return self.get(key)

        if self.level > 0:
            return self.setdefault(key, RecursiveDefaultDict(self.level - 1))

        return self.setdefault(key, dict())


def sliceFromString(string: str, shift: int = 0):
    """Generate a slice from a string, which can represent a slice or an index.

    Parameters
    ----------
    string
        The string representing the slice/index.
    shift
        Potential shift.

    Returns
    -------
    type
        The slice!
    """

    if ":" in string:
        a, b = string.split(":")
        return slice(int(a) + shift, int(b) + shift)
    else:
        return slice(int(string) + shift, int(string) + 1 + shift)


def makeExtractionFunction(expression, symbol="x"):
    """make a simple f(x) expression from string"""
    return lambda x: eval(expression, globals(), {symbol: x})
