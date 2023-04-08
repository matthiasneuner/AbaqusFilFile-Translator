from collections import defaultdict


def RecursiveDefaultDict():
    return defaultdict(RecursiveDefaultDict)


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
