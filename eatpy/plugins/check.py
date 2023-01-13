import datetime

import numpy as np

from ..shared import Plugin, Filter


class Finite(Plugin):
    def before_analysis(
        self,
        time: datetime.datetime,
        state: np.ndarray,
        iobs: np.ndarray,
        obs: np.ndarray,
        obs_sds: np.ndarray,
        filter: Filter,
    ):
        finite = True
        if not np.isfinite(filter.model_states).all():
            self.logger.error("Non-finite values in ensemble state state")
            finite = False
        if not np.isfinite(obs).all():
            self.logger.error("Non-finite values in observations")
            finite = False
        if not np.isfinite(obs_sds).all():
            self.logger.error("Non-finite values in observation standard deviations")
            finite = False
        if not finite:
            raise Exception("non-finite values in model state or observations")
