# the latest flexs (0.2.8) is not compatible with pandas >2.0 as DataFrame.append was depreciated.


import pandas as pd
if not hasattr(pd.DataFrame, "append"):
    def _append_compat(self, other, ignore_index=False, sort=False, **kwargs):
        return pd.concat([self, other], ignore_index=ignore_index, sort=sort)
    pd.DataFrame.append = _append_compat  # type: ignore[attr-defined]