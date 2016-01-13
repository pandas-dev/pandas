
    def _init_mgr(self, mgr, axes=None, dtype=None, copy=False):
        """ passed a manager and a axes dict """
        for a, axe in axes.items():
            if axe is not None:
                mgr = mgr.reindex_axis(
                    axe, axis=self._get_block_manager_axis(a), copy=False)

        # make a copy if explicitly requested
        if copy:
            mgr = mgr.copy()
        if dtype is not None:
            # avoid further copies if we can
            if len(mgr.blocks) > 1 or mgr.blocks[0].values.dtype != dtype:
                mgr = mgr.astype(dtype=dtype)
        return mgr

mgr = self._init_mgr(data, axes=dict(index=index, columns=columns),
                                 dtype=dtype, copy=copy)
