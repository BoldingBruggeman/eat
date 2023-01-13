from typing import Any, MutableMapping
import collections

import numpy as np

from .. import shared


class Spread(shared.Plugin):
    def initialize(self, variables: MutableMapping[str, Any], *args, **kwargs):
        self.variables = collections.OrderedDict(variables)

    def before_analysis(self, *args, **kwargs):
        self.logger.info("Ensemble spread in forecast (root-mean-square differences):")
        for name, info in self.variables.items():
            d = info.get("model_data")
            if d is not None:
                rms = np.sqrt(np.mean(np.var(d, axis=0)))
                self.logger.info("  %s: %.3g %s" % (name, rms, info["units"]))
