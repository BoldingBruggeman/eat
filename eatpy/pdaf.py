from typing import Optional, Sequence, Union
import enum
import os

import numpy as np
from mpi4py import MPI

from . import shared


def import_pdaf():
    global _eat_filter_pdaf
    try:
        from . import _eat_filter_pdaf
    except ImportError as e:
        raise Exception(
            "Import of _eat_filter_pdaf failed."
            " Are MPI, BLAS and LAPACK libraries on your PATH? %s" % e
        )


try:
    import_pdaf()
except Exception:
    if "MKLROOT" in os.environ and hasattr(os, "add_dll_directory"):
        mkldir = os.path.join(os.environ["MKLROOT"], "redist/intel64")
        with os.add_dll_directory(mkldir):
            import_pdaf()
    else:
        raise


class FilterType(enum.IntEnum):
    SEIK = 1
    EnKF = 2
    LSEIK = 3
    ETKF = 4
    LETKF = 5
    ESTKF = 6
    LESTKF = 7
    _3DVar = 13


class CvtHandler(shared.Plugin):
    def __init__(
        self, dim_cvec: Optional[int] = None, dim_cvec_ens: Optional[int] = None
    ):
        self.dim_cvec = dim_cvec
        self.dim_cvec_ens = dim_cvec_ens

    def cvt(self, iter: int, state: np.ndarray, v: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt called but not implemented; state shape = %s, v shape = %s"
            % (state.shape, v.shape,)
        )

    def cvt_adj(self, iter: int, state: np.ndarray, Vv: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_adj called but not implemented; state shape = %s, Vv shape = %s"
            % (state.shape, Vv.shape,)
        )

    def cvt_ens(self, iter: int, state: np.ndarray, v: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_ens called but not implemented; state shape = %s, v shape = %s"
            % (state.shape, v.shape)
        )

    def cvt_adj_ens(self, iter: int, state: np.ndarray, Vv: np.ndarray) -> np.ndarray:
        raise Exception(
            "cvt_adj_ens called but not implemented; state shape = %s, Vv shape = %s"
            % (state.shape, Vv.shape)
        )


class PDAF(shared.Filter):
    """Filter class that wraps [PDAF](https://pdaf.awi.de)."""

    def __init__(
        self,
        filtertype: Union[int, FilterType],
        subtype: int = 0,
        *,
        incremental: int = 0,
        forget: float = 1.0,
        type_forget: int = 0,
        type_trans: int = 0,
        type_sqrt: int = 0,
        type_opt: int = 0,
        rank_analysis_enkf: int = 0,
        beta_3dvar: float = 0.5,
        screen: int = 0,
    ):
        """Initialize PDAF. To get a description of all available option for a
        given filtertype, use subtype = -1
        
        Args:
            filtertype: type of data assimilation filter
            subtype: subtype of data assimilation filter
            incremental: whether to perform incremental analysis
            forget: forgetting factor (usually >0 and <=1)
            type_forget: type of forgetting factor (0: fixed, 1: adaptive)
            type_trans: type of ensemble transformation matrix
                (0: deterministic transformation, 2: use product of 0 with
                random orthonomal matrix with eigenvector (1,...,1)^T)
            type_sqrt: type of transformation matrix square root
                (0: symmetric, 1: Cholesky decomposition)
            type_opt: optimization method (solver) for 3D-Var
                (1: LBFGS, 2: CG+, 3: plain CG)
            rank_analysis_enkf: maximum rank for inversion of HPH^T
                (if 0, HPH is inverted by solving the representer equation)
                (if set to >=ensemble size, it is reset to ensemble size - 1)
            beta_3dvar: weight beta for hybrid 3D-Var
        """
        self.filtertype = FilterType(filtertype)
        self.subtype = subtype
        self.incremental = incremental
        self.forget = forget
        self.type_forget = type_forget
        self.type_trans = type_trans
        self.type_sqrt = type_sqrt
        self.type_opt = type_opt
        self.rank_analysis_enkf = rank_analysis_enkf
        self.beta_3dvar = beta_3dvar
        self.screen = screen

    def initialize(
        self,
        comm: MPI.Comm,
        state_size: int,
        ensemble_size: int,
        plugins: Sequence[shared.Plugin],
    ):
        # Determine if we have a plugin implemnting the CvtHandler API
        # (needed for 3D-Var)
        cvt_handler = None
        for plugin in plugins:
            if isinstance(plugin, CvtHandler):
                cvt_handler = plugin
        if self.filtertype != FilterType._3DVar and cvt_handler is not None:
            raise Exception(
                "One of your plugins implements the CvtHandler routines,"
                " but these are only used when using filtertype 3D-Var (%i)."
                " You are currently using filtertype %s (%i)"
                % (FilterType._3DVar.value, self.filtertype.name, self.filtertype.value)
            )

        if self.filtertype == FilterType.EnKF:
            filter_param_i = [
                state_size,
                ensemble_size,
                self.rank_analysis_enkf,
                self.incremental,
                0,
            ]
            filter_param_r = [self.forget]
        elif self.filtertype == FilterType._3DVar:
            if cvt_handler is None:
                raise Exception(
                    "To use filtertype 3D-Var (%i),"
                    " one of your plugins must derive from eatpy.pdaf.CvtHandler"
                    % FilterType._3DVar.value
                )
            _eat_filter_pdaf.cvt_handler = cvt_handler

            if self.subtype in (0, 6, 7) and cvt_handler.dim_cvec is None:
                raise Exception(
                    "For parameterized or hybrid 3D-Var, the attribute dim_cvec"
                    " must be set in your CvtHandler plugin."
                )
            filter_param_i = [
                state_size,
                ensemble_size,
                self.type_opt,
                cvt_handler.dim_cvec or -1,
                cvt_handler.dim_cvec_ens or ensemble_size,
            ]
            filter_param_r = [self.forget, self.beta_3dvar]
        else:
            filter_param_i = [
                state_size,
                ensemble_size,
                0,
                self.incremental,
                self.type_forget,
                self.type_trans,
                self.type_sqrt,
            ]
            filter_param_r = [self.forget]
        filter_param_i = np.array(filter_param_i, dtype=np.intc)
        filter_param_r = np.array(filter_param_r, dtype=np.double)
        # self.model_states = _eat_filter_pdaf.initialize(comm, state_size, ensemble_size)
        self.model_states = _eat_filter_pdaf.initialize_with_params(
            self.filtertype,
            self.subtype,
            filter_param_i,
            filter_param_r,
            comm,
            self.screen,
        )

    def assimilate(self, iobs: np.ndarray, obs: np.ndarray, sds: np.ndarray):
        iobs = iobs + 1  # convert to 1-based indices for Fortran/PDAF
        _eat_filter_pdaf.assimilate(iobs, obs, sds)

    def finalize(self):
        _eat_filter_pdaf.finalize()

