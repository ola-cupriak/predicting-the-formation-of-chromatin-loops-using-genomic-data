from kedro.io import AbstractDataSet


class PathDataSet(AbstractDataSet):
    """
    ``PathDataSet`` loads data from a given filepath as path.
    """

    def __init__(self, filepath: str):
        """Creates a new instance of bigWigDataSet to load / save image data for given filepath.

        Args:
            filepath: The location of the bigWig file to load / save data.
        """
        self._filepath = filepath

    def _load(self):
        """
        Return the filepath.
        """
        return self._filepath

    def _save(self) -> None:
        return None

    def _describe(self) -> None:
        """Returns a dict that describes the attributes of the dataset."""
        return None