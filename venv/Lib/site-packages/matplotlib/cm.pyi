from collections.abc import Iterator, Mapping
from matplotlib import colors
from matplotlib.colorizer import _ScalarMappable


class ColormapRegistry(Mapping[str, colors.Colormap]):
    def __init__(self, cmaps: Mapping[str, colors.Colormap]) -> None: ...
    def __getitem__(self, item: str) -> colors.Colormap: ...
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...
    def __call__(self) -> list[str]: ...
    def register(
        self, cmap: colors.Colormap, *, name: str | None = ..., force: bool = ...
    ) -> None: ...
    def unregister(self, name: str) -> None: ...
    def get_cmap(self, cmap: str | colors.Colormap) -> colors.Colormap: ...

_colormaps: ColormapRegistry = ...
_multivar_colormaps: ColormapRegistry = ...
_bivar_colormaps: ColormapRegistry = ...

ScalarMappable = _ScalarMappable

magma: colors.Colormap
inferno: colors.Colormap
plasma: colors.Colormap
viridis: colors.Colormap
cividis: colors.Colormap
twilight: colors.Colormap
twilight_shifted: colors.Colormap
turbo: colors.Colormap
berlin: colors.Colormap
managua: colors.Colormap
vanimo: colors.Colormap
Blues: colors.Colormap
BrBG: colors.Colormap
BuGn: colors.Colormap
BuPu: colors.Colormap
CMRmap: colors.Colormap
GnBu: colors.Colormap
Greens: colors.Colormap
Greys: colors.Colormap
OrRd: colors.Colormap
Oranges: colors.Colormap
PRGn: colors.Colormap
PiYG: colors.Colormap
PuBu: colors.Colormap
PuBuGn: colors.Colormap
PuOr: colors.Colormap
PuRd: colors.Colormap
Purples: colors.Colormap
RdBu: colors.Colormap
RdGy: colors.Colormap
RdPu: colors.Colormap
RdYlBu: colors.Colormap
RdYlGn: colors.Colormap
Reds: colors.Colormap
Spectral: colors.Colormap
Wistia: colors.Colormap
YlGn: colors.Colormap
YlGnBu: colors.Colormap
YlOrBr: colors.Colormap
YlOrRd: colors.Colormap
afmhot: colors.Colormap
autumn: colors.Colormap
binary: colors.Colormap
bone: colors.Colormap
brg: colors.Colormap
bwr: colors.Colormap
cool: colors.Colormap
coolwarm: colors.Colormap
copper: colors.Colormap
cubehelix: colors.Colormap
flag: colors.Colormap
gist_earth: colors.Colormap
gist_gray: colors.Colormap
gist_heat: colors.Colormap
gist_ncar: colors.Colormap
gist_rainbow: colors.Colormap
gist_stern: colors.Colormap
gist_yarg: colors.Colormap
gnuplot: colors.Colormap
gnuplot2: colors.Colormap
gray: colors.Colormap
hot: colors.Colormap
hsv: colors.Colormap
jet: colors.Colormap
nipy_spectral: colors.Colormap
ocean: colors.Colormap
pink: colors.Colormap
prism: colors.Colormap
rainbow: colors.Colormap
seismic: colors.Colormap
spring: colors.Colormap
summer: colors.Colormap
terrain: colors.Colormap
winter: colors.Colormap
Accent: colors.Colormap
Dark2: colors.Colormap
okabe_ito: colors.Colormap
Paired: colors.Colormap
Pastel1: colors.Colormap
Pastel2: colors.Colormap
Set1: colors.Colormap
Set2: colors.Colormap
Set3: colors.Colormap
tab10: colors.Colormap
tab20: colors.Colormap
tab20b: colors.Colormap
tab20c: colors.Colormap
grey: colors.Colormap
gist_grey: colors.Colormap
gist_yerg: colors.Colormap
Grays: colors.Colormap
# Reversed colormaps
magma_r: colors.Colormap
inferno_r: colors.Colormap
plasma_r: colors.Colormap
viridis_r: colors.Colormap
cividis_r: colors.Colormap
twilight_r: colors.Colormap
twilight_shifted_r: colors.Colormap
turbo_r: colors.Colormap
berlin_r: colors.Colormap
managua_r: colors.Colormap
vanimo_r: colors.Colormap
Blues_r: colors.Colormap
BrBG_r: colors.Colormap
BuGn_r: colors.Colormap
BuPu_r: colors.Colormap
CMRmap_r: colors.Colormap
GnBu_r: colors.Colormap
Greens_r: colors.Colormap
Greys_r: colors.Colormap
OrRd_r: colors.Colormap
Oranges_r: colors.Colormap
PRGn_r: colors.Colormap
PiYG_r: colors.Colormap
PuBu_r: colors.Colormap
PuBuGn_r: colors.Colormap
PuOr_r: colors.Colormap
PuRd_r: colors.Colormap
Purples_r: colors.Colormap
RdBu_r: colors.Colormap
RdGy_r: colors.Colormap
RdPu_r: colors.Colormap
RdYlBu_r: colors.Colormap
RdYlGn_r: colors.Colormap
Reds_r: colors.Colormap
Spectral_r: colors.Colormap
Wistia_r: colors.Colormap
YlGn_r: colors.Colormap
YlGnBu_r: colors.Colormap
YlOrBr_r: colors.Colormap
YlOrRd_r: colors.Colormap
afmhot_r: colors.Colormap
autumn_r: colors.Colormap
binary_r: colors.Colormap
bone_r: colors.Colormap
brg_r: colors.Colormap
bwr_r: colors.Colormap
cool_r: colors.Colormap
coolwarm_r: colors.Colormap
copper_r: colors.Colormap
cubehelix_r: colors.Colormap
flag_r: colors.Colormap
gist_earth_r: colors.Colormap
gist_gray_r: colors.Colormap
gist_heat_r: colors.Colormap
gist_ncar_r: colors.Colormap
gist_rainbow_r: colors.Colormap
gist_stern_r: colors.Colormap
gist_yarg_r: colors.Colormap
gnuplot_r: colors.Colormap
gnuplot2_r: colors.Colormap
gray_r: colors.Colormap
hot_r: colors.Colormap
hsv_r: colors.Colormap
jet_r: colors.Colormap
nipy_spectral_r: colors.Colormap
ocean_r: colors.Colormap
pink_r: colors.Colormap
prism_r: colors.Colormap
rainbow_r: colors.Colormap
seismic_r: colors.Colormap
spring_r: colors.Colormap
summer_r: colors.Colormap
terrain_r: colors.Colormap
winter_r: colors.Colormap
Accent_r: colors.Colormap
Dark2_r: colors.Colormap
okabe_ito_r: colors.Colormap
Paired_r: colors.Colormap
Pastel1_r: colors.Colormap
Pastel2_r: colors.Colormap
Set1_r: colors.Colormap
Set2_r: colors.Colormap
Set3_r: colors.Colormap
tab10_r: colors.Colormap
tab20_r: colors.Colormap
tab20b_r: colors.Colormap
tab20c_r: colors.Colormap
grey_r: colors.Colormap
gist_grey_r: colors.Colormap
gist_yerg_r: colors.Colormap
Grays_r: colors.Colormap
