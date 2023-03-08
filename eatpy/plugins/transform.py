from typing import Any, List, MutableMapping
import datetime

import numpy as np

from .. import shared


class Log(shared.Plugin):
    def __init__(
        self,
        *variable_names,
        transform_obs: bool = True,
        minimum: float = -np.inf,
        log10: bool = True
    ):
        self.variable_names = frozenset(variable_names)
        self.variable_metadata: List[Any] = []
        self.transform_obs = transform_obs
        self.minimum = minimum
        self.log10 = log10
        self.forward = np.log10 if self.log10 else np.log
        self.backward = (lambda x: 10.0 ** x) if self.log10 else np.exp

    def initialize(self, variables: MutableMapping[str, Any], *args, **kwargs):
        for name in self.variable_names:
            self.variable_metadata.append(variables[name])

    def before_analysis(
        self,
        time: datetime.datetime,
        state: np.ndarray,
        iobs: np.ndarray,
        obs: np.ndarray,
        obs_sds: np.ndarray,
        *args,
        **kwargs
    ):
        for metadata in self.variable_metadata:
            affected_obs = (iobs >= metadata["start"]) & (iobs < metadata["stop"])
            if self.transform_obs and affected_obs.any():
                # Transform mean and sd of observations,
                # assuming their distribution is log-normal
                mean = np.maximum(obs[affected_obs], self.minimum)
                sd = obs_sds[affected_obs]
                sigma2 = np.log((sd / mean) ** 2 + 1.0)
                mu = np.log(mean) - 0.5 * sigma2
                sigma = np.sqrt(sigma2)
                if self.log10:
                    mu /= np.log(10.0)
                    sigma /= np.log(10.0)
                obs[affected_obs] = mu
                obs_sds[affected_obs] = sigma
            metadata["data"][...] = self.forward(
                np.maximum(metadata["data"], self.minimum)
            )

    def after_analysis(self, *args, **kwargs):
        for metadata in self.variable_metadata:
            metadata["data"][...] = self.backward(metadata["data"])
