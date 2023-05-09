from dask.distributed import Client


class DaskClient(Client):
    def __init__(self):
        super().__init__()

    def cleanup(self):
        self.close()
