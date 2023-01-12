from typing import MutableMapping, Any

from .select import Select
from . import transform
from . import output
from ..shared import Plugin


class TestPlugin(Plugin):
    def initialize(self, variables: MutableMapping[str, Any], ensemble_size: int):
        self.logger.info("initialize")

    def before_analysis(self, *args, **kwargs):
        self.logger.info("before_analysis")

    def after_analysis(self, *args, **kwargs):
        self.logger.info("after_analysis")

    def finalize(self):
        self.logger.info("finalize")
