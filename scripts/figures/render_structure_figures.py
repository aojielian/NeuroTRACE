from pathlib import Path
import os,importlib.util
from matplotlib.text import Text
from PIL import Image
root=Path(os.environ.get('NEUROTRACE_PROJECT_ROOT','.')).resolve()
output=Path(os.environ.get('NEUROTRACE_FIGURE_OUTPUT',str(root/'generated_figures')));output.mkdir(parents=True,exist_ok=True)
s=importlib.util.spec_from_file_location('structures',root/'scripts/figures/current_structures.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def export(f,name):
 newname=name.replace('Figure_S','FigureS');scale=180/(8.3*25.4);w,h=f.get_size_inches();f.set_size_inches(w*scale,h*scale)
 for t in f.findobj(Text):t.set_fontsize(t.get_fontsize()*scale)
 out=output/newname;f.savefig(out.with_suffix('.pdf'),metadata={'Title':newname,'Creator':'Scientific plotting'});f.savefig(out.with_suffix('.png'),dpi=600)
 with Image.open(out.with_suffix('.png')) as im:im.convert('RGB').save(out.with_suffix('.tiff'),dpi=(600,600),compression='tiff_lzw')
 m.plt.close(f);print(newname,flush=True)
m.export=export
for n in ['s1','s4','s6']:getattr(m,n)()
