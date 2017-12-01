---
layout: post
title: Interactive Genome plots using bokeh
---



{% highlight python %}
from bokeh.plotting import figure, output_notebook, show, gridplot
from bokeh.models.tickers import FixedTicker
from bokeh.models import NumeralTickFormatter, Range1d, LinearAxis, ColumnDataSource, HoverTool
{% endhighlight %}


{% highlight python %}
import sys
import os
import pyfasta
import pandas as pd
import numpy as np
import allel
import functools
{% endhighlight %}


{% highlight python %}
# output to static HTML file
output_notebook()
{% endhighlight %}



    <div class="bk-root">
        <a href="https://bokeh.pydata.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
        <span id="c175191b-8022-49b5-b522-81e00aa42e21">Loading BokehJS ...</span>
    </div>




### Load ag1k data


{% highlight python %}
sys.path.insert(0, '../../../selection_paper/agam-report-base/src/python')
ag1k_dir = '/kwiat/vector/ag1000g/release'
from ag1k import phase1_ar3, phase1_ar31
{% endhighlight %}


{% highlight python %}
phase1_ar3.init(os.path.join(ag1k_dir, 'phase1.AR3'))
{% endhighlight %}


{% highlight python %}
chromosomes = "2R", "2L", "3R", "3L", "X"
{% endhighlight %}

### Define genome plot class


{% highlight python %}
class GenomePlot:
    
    genome = None
    contigs = None
    tools = "pan,wheel_zoom,box_zoom,save,reset"
    layout = None
    pfunc = None
    min_border = 50
    plot_width_per_mb = 6
    
    def __init__(self, fasta, contigs=None, layout=None, share_y=True, pfunc=None):
        
        self.genome = pyfasta.Fasta(fasta)
        
        if contigs is None:
            self.contigs = list(self.genome.keys())
        else:
            self.contigs = contigs
            
        # handle layout
        if layout is not None:
            self.layout = self.parse_layout(layout)
        else:
            self.layout = self.auto_layout()
        
        if pfunc is not None:
            self.pfunc = pfunc
                
    def parse_layout(self, lstring):
        layout = list()
        assert type(lstring) == str, "layout must be a string"
        q = list(chromosomes)
        q.reverse()
            
        rows = lstring.split("|")
        for row in rows:
            r = list()
            for x in row:
                assert x in "ox.", "Only chars o, x, . are allowed"
                if x == ".":
                    r.append(None)
                else:
                    r.append(q.pop())
            layout.append(r)
        return layout
    
    def auto_layout():
        pass

    def apply(self, func, **kwargs):
        
        # create a figure with specified layout        
        d = [[] for i in range(len(self.layout))]
        
        for i, row in enumerate(self.layout):
            
            for chrom in row:
                
                if chrom is not None:
                    
                    csize = len(self.genome[chrom])
                    px = int(csize * 1e-6 * self.plot_width_per_mb) + self.min_border
                    
                    try:
                        yrange = s1.y_range
                    except NameError:
                        yrange = None

                    s1 = figure(width=px, 
                                min_border_left=self.min_border,
                                plot_height=250, 
                                tools=self.tools,
                                title=chrom,
                                y_range=yrange,
                                x_range=(1, csize))

                    s1.xaxis.ticker = FixedTicker(ticks=np.arange(0, csize, 1e7))
                    s1.xaxis[0].formatter = NumeralTickFormatter(format="0a.0")
                    
                    # handle general plot things specific to genome not data
                    if self.pfunc is not None:
                        self.pfunc(chrom, s1)

                    # function lives here
                    func(chrom, s1, **kwargs)
                    d[i].append(s1)

                else:
                    d[i].append(None)
            
            
        # NEW: put the subplots in a gridplot
        p = gridplot(d, toolbar_location="left", sizing_mode='fixed', plot_width=None)
        show(p)
{% endhighlight %}


{% highlight python %}
# Custom changes to axes that pertain to genome not data
# eg in this case move left chromosome arm axes to RHS
def custom_plot(chrom, subplot):
    if chrom.endswith("L"):
        subplot.yaxis.visible = False
        subplot.add_layout(LinearAxis(), 'right')
{% endhighlight %}

### Calculate $\theta$, $\Pi$, and Tajima's D


{% highlight python %}
@functools.lru_cache()
def calculate_summary_stats(chrom, pop, window_size=100000):
    
    ix = phase1_ar3.df_samples.query("population == @pop").index
    accessibility = phase1_ar3.accessibility[chrom]["is_accessible"][:]
    
    pos = allel.SortedIndex(phase1_ar3.callset_pass[chrom]["variants/POS"][:])
    eqw = allel.equally_accessible_windows(accessibility, window_size)
    g = allel.GenotypeChunkedArray(phase1_ar3.callset_pass[chrom]["calldata/genotype"]).take(ix, axis=1)
    ac = g.count_alleles()
    
    theta, wins, nb, counts = allel.stats.windowed_watterson_theta(pos, ac, windows=eqw, is_accessible=accessibility)
    
    pi, wins, nb, counts = allel.stats.windowed_diversity(pos, ac, windows=eqw, is_accessible=accessibility)
    
    tajD, wins, counts = allel.stats.windowed_tajima_d(pos, ac, windows=eqw)
    
    df = pd.DataFrame.from_dict({"start": eqw[:, 0], 
                                 "stop": eqw[:, 1], 
                                 "diversity": pi, 
                                 "tajimaD": tajD, 
                                 "theta": theta})

    df["midpoint"] = eqw.mean(1)
    
    return df
{% endhighlight %}


{% highlight python %}
stats = {c: calculate_summary_stats(chrom=c, pop="BFS", window_size=100000) for c in chromosomes}
{% endhighlight %}

### Use a .gff3 file to annotate above windows


{% highlight python %}
gff3 = allel.FeatureTable.from_gff3(phase1_ar3.geneset_agamp42_fn, attributes=["ID"])
{% endhighlight %}


{% highlight python %}
gff3
{% endhighlight %}




<div class="allel allel-DisplayAsTable"><span>&lt;FeatureTable shape=(175804,) dtype=(numpy.record, [('seqid', 'O'), ('source', 'O'), ('type', 'O'), ('start', '&lt;i8'), ('end', '&lt;i8'), ('score', '&lt;f8'), ('strand', 'O'), ('phase', '&lt;i8'), ('ID', 'O')])&gt;</span><table><thead><tr><th></th><th style="text-align: center">seqid</th><th style="text-align: center">source</th><th style="text-align: center">type</th><th style="text-align: center">start</th><th style="text-align: center">end</th><th style="text-align: center">score</th><th style="text-align: center">strand</th><th style="text-align: center">phase</th><th style="text-align: center">ID</th></tr></thead><tbody><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">0</th><td style="text-align: center">2L</td><td style="text-align: center">VectorBase</td><td style="text-align: center">contig</td><td style="text-align: center">1</td><td style="text-align: center">49364325</td><td style="text-align: center">-1.0</td><td style="text-align: center">.</td><td style="text-align: center">-1</td><td style="text-align: center">2L</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">1</th><td style="text-align: center">2L</td><td style="text-align: center">VectorBase</td><td style="text-align: center">exon</td><td style="text-align: center">157348</td><td style="text-align: center">157623</td><td style="text-align: center">-1.0</td><td style="text-align: center">-</td><td style="text-align: center">-1</td><td style="text-align: center">AGAP004677-RB-E4A</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">2</th><td style="text-align: center">2L</td><td style="text-align: center">VectorBase</td><td style="text-align: center">exon</td><td style="text-align: center">157348</td><td style="text-align: center">157623</td><td style="text-align: center">-1.0</td><td style="text-align: center">-</td><td style="text-align: center">-1</td><td style="text-align: center">AGAP004677-RB-E4B</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">...</th><td style="text-align: center" colspan="10">...</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">175801</th><td style="text-align: center">X</td><td style="text-align: center">VectorBase</td><td style="text-align: center">gene</td><td style="text-align: center">24338771</td><td style="text-align: center">24340371</td><td style="text-align: center">-1.0</td><td style="text-align: center">+</td><td style="text-align: center">-1</td><td style="text-align: center">AGAP013609</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">175802</th><td style="text-align: center">X</td><td style="text-align: center">VectorBase</td><td style="text-align: center">rRNA</td><td style="text-align: center">24338771</td><td style="text-align: center">24340371</td><td style="text-align: center">-1.0</td><td style="text-align: center">+</td><td style="text-align: center">-1</td><td style="text-align: center">AGAP013609-RA</td></tr><tr><th style="text-align: center; background-color: white; border-right: 1px solid black; ">175803</th><td style="text-align: center">Y_unplaced</td><td style="text-align: center">VectorBase</td><td style="text-align: center">contig</td><td style="text-align: center">1</td><td style="text-align: center">237045</td><td style="text-align: center">-1.0</td><td style="text-align: center">.</td><td style="text-align: center">-1</td><td style="text-align: center">Y_unplaced</td></tr></tbody></table></div>




{% highlight python %}
annotated_data = {}

# annotate these data
for chrom in chromosomes:

    d = stats[chrom].copy()
    
    # extract the relevant seq id and use pandas interval indexing
    features = pd.DataFrame(gff3.query("seqid == '{0}'".format(chrom)).values)
    features.index = pd.IntervalIndex.from_arrays(features.start, features.end, closed="both")

    # logic to extract relevant rows, filter by annot type, drop duplicates and join ID column
    # it would be slightly more efficient to do these both in a single call, but it's fast/readable so we'll let it slide
    d["gene"] = d.apply(
        lambda y: ", ".join(features.loc[[y.start, y.stop]].query("type == 'gene'").ID.drop_duplicates()), 1)

    annotated_data[chrom] = d
{% endhighlight %}


{% highlight python %}
annotated_data["X"].head()
{% endhighlight %}




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>diversity</th>
      <th>start</th>
      <th>stop</th>
      <th>tajimaD</th>
      <th>theta</th>
      <th>midpoint</th>
      <th>gene</th>
      <th>contig</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.001781</td>
      <td>25</td>
      <td>132324</td>
      <td>-2.190824</td>
      <td>0.005357</td>
      <td>66174.5</td>
      <td>AGAP000011</td>
      <td>X</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.003892</td>
      <td>132325</td>
      <td>246994</td>
      <td>-2.041654</td>
      <td>0.010288</td>
      <td>189659.5</td>
      <td>AGAP000011, AGAP000018</td>
      <td>X</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.005919</td>
      <td>246995</td>
      <td>368599</td>
      <td>-2.176849</td>
      <td>0.017548</td>
      <td>307797.0</td>
      <td>AGAP000018</td>
      <td>X</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.006699</td>
      <td>368600</td>
      <td>487739</td>
      <td>-2.297078</td>
      <td>0.022272</td>
      <td>428169.5</td>
      <td></td>
      <td>X</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.009483</td>
      <td>487740</td>
      <td>625447</td>
      <td>-2.252033</td>
      <td>0.030146</td>
      <td>556593.5</td>
      <td></td>
      <td>X</td>
    </tr>
  </tbody>
</table>
</div>




{% highlight python %}
def stat_plot(chrom, subplot, stat="diversity", label=None):
    
    source = ColumnDataSource(annotated_data[chrom])
    
    if label is None:
        label = stat
    
    hover = HoverTool(
        tooltips=[
            ("Position", "@start{0a.000}-@stop{0a.000}"),
            (label, "$y"),
            ("contig", chrom), 
            ("gene(s)", "@gene")],
        mode="mouse")
    
    subplot.add_tools(hover)
    
    subplot.circle("midpoint", 
                   stat,
                   source=source,
                   size=3, 
                   color="navy", 
                   alpha=0.5)

    subplot.yaxis[0].axis_label = label
    subplot.xaxis[0].axis_label = "Genomic Position (bp)"
{% endhighlight %}


{% highlight python %}
gf = GenomePlot(fasta=phase1_ar3.genome_fn, contigs=("2R", "2L", "3R", "3L", "X"), layout="oo|ooo", pfunc=custom_plot)
{% endhighlight %}


{% highlight python %}
gf.apply(stat_plot, stat="tajimaD", label="Tajima'sD")
{% endhighlight %}



<div class="bk-root">
    <div class="bk-plotdiv" id="09cf6405-8ee6-4dfb-a6e6-9af47d39ddbe"></div>
</div>





{% highlight python %}
gf.apply(stat_plot, stat="theta", label="\u0398W")
{% endhighlight %}



<div class="bk-root">
    <div class="bk-plotdiv" id="d9fb87f4-9015-4b00-ae83-db8cb22f5af4"></div>
</div>





{% highlight python %}
gf.apply(stat_plot, stat="diversity", label="\u03A0")
{% endhighlight %}



<div class="bk-root">
    <div class="bk-plotdiv" id="23680b8a-f988-4f14-b7ab-53c58c923b1d"></div>
</div>





{% highlight python %}

{% endhighlight %}
