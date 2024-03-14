from tools.data_object import DataObject

from plot_inputs import PlotInputData, SLICE_TYPE


class SingleLineInputData(PlotInputData):
    def __init__(self, data_object: DataObject, coast_index: int,
                 indicator_index: int):
        super().__init__(data_object)
        self._coast_index: int = coast_index
        self._indicator_index: int = indicator_index

    def get_time_steps(self):
        return self._data_object.get_time_steps()

    def get_values(self, time_steps: SLICE_TYPE, basin: SLICE_TYPE,
                   level_index: SLICE_TYPE):
        return self._data_object.get_values(
            time_steps,
            basin,
            level_index,
            coasts=self._coast_index,
            indicator=self._indicator_index
        )

    def __enter__(self):
        self._data_object.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._data_object.__exit__(exc_type, exc_val, exc_tb)
        return
