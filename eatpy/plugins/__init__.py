from typing import MutableMapping, Any

from ..shared import Plugin
from . import transform
from . import output
from . import select
from . import report
from . import check
from .select import Select


class TestPlugin(Plugin):
    def initialize(self, variables: MutableMapping[str, Any], ensemble_size: int):
        self.logger.info(f"initialize (ensemble size {ensemble_size})")

    def before_analysis(self, *args, **kwargs):
        self.logger.info("before_analysis")

    def after_analysis(self, *args, **kwargs):
        self.logger.info("after_analysis")

    def finalize(self):
        self.logger.info("finalize")
