
<html>
<h1>Moving Average Plot R Code</h2>
<p>Created by: <b>Awni Mousa</b><p>
<p><i>July 29th, 2015</i></p>
</head>

<body>

<p>The following R function will apply a sliding window of moving average method given a window size and steps.</p>
<p>As input, a two column <b>data frame</b>, <i>df</i>, of a response and terms, is used, where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. Example, response = log2FC, predictor = geneLength. Optional arguments can be changed to control the window size and the step value.
<p>The <i>winSize</i>, an integer > 1, controls the desired window size; default 200.</p>
<p>The <i>step</i> controls the distance between windows and is an integer <= <i>winSize</i>; default 40.</p>
<p>The output is a <i>LIST</i> of 2 elements:
<ul>
<li>a 4 column data frame:
<ul style="list-style-type:circle">
<li>predMean - the mean of the predictor terms within a window</li>
<li>respAvg - the mean of the response within a window (the moving average)</li>
<li>ucl, lcl - the 95% CI upper and lower levels</li></ul>
<li>a plotting object</li></ul>
<p>These values can be assigned to a variable to be used in the environment by calling the first element in the <i>list[[1]]</i> to get the data frame with the values, and by calling the second element in the <i>list[[2]]</i> to initiate the plot (<i>see example below</i>).</p>

<h2>R code</h2>

<div class="chunk" id="unnamed-chunk-1"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl std">movingAverage.AM</span> <span class="hl kwb">&lt;-</span> <span class="hl kwa">function</span><span class="hl std">(</span> <span class="hl kwc">df</span><span class="hl std">,</span> <span class="hl kwc">winSize</span><span class="hl std">=</span><span class="hl num">200</span><span class="hl std">,</span> <span class="hl kwc">step</span><span class="hl std">=</span><span class="hl num">40</span><span class="hl std">){</span>
  <span class="hl kwd">colnames</span><span class="hl std">(df)</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">c</span><span class="hl std">(</span><span class="hl str">&quot;response&quot;</span><span class="hl std">,</span><span class="hl str">&quot;predictor&quot;</span><span class="hl std">)</span>
  <span class="hl std">meanValue</span> <span class="hl kwb">&lt;-</span> <span class="hl kwa">NULL</span>
  <span class="hl std">l</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">seq</span><span class="hl std">(</span><span class="hl num">1</span><span class="hl std">,</span><span class="hl kwd">length</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response),step)</span>
  <span class="hl kwa">for</span><span class="hl std">(i</span> <span class="hl kwa">in</span> <span class="hl std">l){</span>
    <span class="hl kwa">if</span><span class="hl std">(i</span><span class="hl opt">+</span><span class="hl std">winSize</span><span class="hl opt">-</span><span class="hl num">1</span> <span class="hl opt">&gt;</span> <span class="hl kwd">length</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response)){</span>
      <span class="hl std">i1</span> <span class="hl kwb">=</span> <span class="hl kwd">length</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response) }</span>
    <span class="hl kwa">else</span> <span class="hl std">i1</span> <span class="hl kwb">=</span> <span class="hl std">i</span><span class="hl opt">+</span><span class="hl std">winSize</span><span class="hl opt">-</span><span class="hl num">1</span>
    <span class="hl com"># mean</span>
    <span class="hl std">respAvg</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">mean</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response[i</span><span class="hl opt">:</span><span class="hl std">i1])</span>
    <span class="hl com"># 95% CI</span>
    <span class="hl std">ucl</span> <span class="hl kwb">&lt;-</span> <span class="hl std">respAvg</span> <span class="hl opt">+</span> <span class="hl num">1.96</span> <span class="hl opt">*</span> <span class="hl kwd">sd</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response[i</span><span class="hl opt">:</span><span class="hl std">i1]</span><span class="hl opt">/</span><span class="hl kwd">sqrt</span><span class="hl std">(</span><span class="hl kwd">length</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response[i</span><span class="hl opt">:</span><span class="hl std">i1])))</span>
    <span class="hl std">lcl</span> <span class="hl kwb">&lt;-</span> <span class="hl std">respAvg</span> <span class="hl opt">-</span> <span class="hl num">1.96</span> <span class="hl opt">*</span> <span class="hl kwd">sd</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response[i</span><span class="hl opt">:</span><span class="hl std">i1]</span><span class="hl opt">/</span><span class="hl kwd">sqrt</span><span class="hl std">(</span><span class="hl kwd">length</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">response[i</span><span class="hl opt">:</span><span class="hl std">i1])))</span>
    <span class="hl com"># average gene length in bin</span>
    <span class="hl std">predMean</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">mean</span><span class="hl std">(df</span><span class="hl opt">$</span><span class="hl std">predictor[i</span><span class="hl opt">:</span><span class="hl std">i1])</span>
    <span class="hl std">meanValue.i</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">cbind</span><span class="hl std">(predMean,respAvg,ucl,lcl)</span>
    <span class="hl std">meanValue</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">rbind</span><span class="hl std">(meanValue, meanValue.i)</span>
  <span class="hl std">}</span>
  <span class="hl std">meanValue</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">as.data.frame</span><span class="hl std">(meanValue)</span>
  <span class="hl com"># make plot</span>
  <span class="hl kwd">library</span><span class="hl std">(ggplot2)</span>
  <span class="hl std">p</span> <span class="hl kwb">=</span> <span class="hl kwd">ggplot</span><span class="hl std">(df,</span> <span class="hl kwd">aes</span><span class="hl std">(</span><span class="hl kwc">y</span><span class="hl std">=response,</span> <span class="hl kwc">x</span><span class="hl std">=predictor))</span>
  <span class="hl std">P</span> <span class="hl kwb">&lt;-</span> <span class="hl std">p</span> <span class="hl opt">+</span> <span class="hl kwd">scale_x_log10</span><span class="hl std">(</span><span class="hl kwc">breaks</span><span class="hl std">=</span><span class="hl num">10</span><span class="hl opt">^</span><span class="hl std">(</span><span class="hl kwd">seq</span><span class="hl std">(</span><span class="hl opt">-</span><span class="hl num">1</span><span class="hl std">,</span><span class="hl num">10</span><span class="hl std">)),</span><span class="hl kwc">labels</span><span class="hl std">=</span><span class="hl num">10</span><span class="hl opt">^</span><span class="hl std">(</span><span class="hl kwd">seq</span><span class="hl std">(</span><span class="hl opt">-</span><span class="hl num">1</span><span class="hl std">,</span><span class="hl num">10</span><span class="hl std">)))</span> <span class="hl opt">+</span>
    <span class="hl kwd">geom_point</span><span class="hl std">(</span><span class="hl kwc">color</span><span class="hl std">=</span><span class="hl str">&quot;grey64&quot;</span><span class="hl std">,</span> <span class="hl kwc">alpha</span><span class="hl std">=</span><span class="hl num">1</span><span class="hl opt">/</span><span class="hl num">2</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">geom_abline</span><span class="hl std">(</span><span class="hl kwc">intercept</span> <span class="hl std">=</span> <span class="hl num">0</span><span class="hl std">,</span> <span class="hl kwc">slope</span> <span class="hl std">=</span> <span class="hl num">0</span><span class="hl std">,</span> <span class="hl kwc">lwd</span> <span class="hl std">=</span> <span class="hl num">1</span><span class="hl std">,</span> <span class="hl kwc">col</span> <span class="hl std">=</span> <span class="hl str">&quot;#009E73&quot;</span><span class="hl std">,</span> <span class="hl kwc">alpha</span><span class="hl std">=</span><span class="hl num">1</span><span class="hl opt">/</span><span class="hl num">2</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">annotation_logticks</span><span class="hl std">(</span><span class="hl kwc">sides</span><span class="hl std">=</span><span class="hl str">&quot;b&quot;</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">geom_smooth</span><span class="hl std">(</span><span class="hl kwd">aes</span><span class="hl std">(</span><span class="hl kwc">x</span><span class="hl std">=predMean,</span> <span class="hl kwc">y</span><span class="hl std">=respAvg,</span> <span class="hl kwc">ymin</span><span class="hl std">=lcl,</span> <span class="hl kwc">ymax</span><span class="hl std">=ucl),</span> <span class="hl kwc">data</span><span class="hl std">=meanValue,</span>
                <span class="hl kwc">stat</span><span class="hl std">=</span><span class="hl str">&quot;identity&quot;</span><span class="hl std">,</span> <span class="hl kwc">color</span><span class="hl std">=</span><span class="hl str">&quot;red&quot;</span><span class="hl std">,</span> <span class="hl kwc">fill</span><span class="hl std">=</span><span class="hl str">&quot;red&quot;</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">geom_rug</span><span class="hl std">(</span><span class="hl kwc">position</span><span class="hl std">=</span><span class="hl str">&quot;jitter&quot;</span><span class="hl std">,</span> <span class="hl kwc">size</span><span class="hl std">=</span><span class="hl num">0.05</span><span class="hl std">,</span> <span class="hl kwc">color</span><span class="hl std">=</span><span class="hl str">&quot;magenta&quot;</span><span class="hl std">,</span> <span class="hl kwc">alpha</span><span class="hl std">=</span><span class="hl num">1</span><span class="hl opt">/</span><span class="hl num">8</span><span class="hl std">,</span>
             <span class="hl kwc">sides</span><span class="hl std">=</span><span class="hl str">&quot;tr&quot;</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">xlab</span><span class="hl std">(</span><span class="hl str">&quot;predictor&quot;</span><span class="hl std">)</span> <span class="hl opt">+</span> <span class="hl kwd">ylab</span><span class="hl std">(</span><span class="hl str">&quot;response&quot;</span><span class="hl std">)</span> <span class="hl opt">+</span>
    <span class="hl kwd">ggtitle</span><span class="hl std">(</span><span class="hl kwd">bquote</span><span class="hl std">(</span> <span class="hl kwd">atop</span><span class="hl std">(</span><span class="hl str">&quot;Moving Average&quot;</span><span class="hl std">,</span> <span class="hl kwd">paste</span><span class="hl std">(</span><span class="hl str">&quot;Window Size &quot;</span><span class="hl std">,</span><span class="hl kwd">.</span><span class="hl std">(winSize),</span>
                                                 <span class="hl str">&quot; Step &quot;</span><span class="hl std">,</span> <span class="hl kwd">.</span><span class="hl std">(step) ))))</span>
  <span class="hl kwd">return</span><span class="hl std">(</span><span class="hl kwd">list</span><span class="hl std">(</span><span class="hl kwc">MovingAverage</span><span class="hl std">=meanValue,</span> <span class="hl kwc">plot</span><span class="hl std">=P))</span>
<span class="hl std">}</span>
</pre></div>
</div></div>

<h3>Usage and example plot</h3>

<div class="chunk" id="unnamed-chunk-2"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl std">df</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">read.table</span><span class="hl std">(</span><span class="hl str">&quot;df.txt&quot;</span><span class="hl std">,</span><span class="hl kwc">header</span> <span class="hl std">= T,</span> <span class="hl kwc">sep</span> <span class="hl std">=</span> <span class="hl str">&quot;\t&quot;</span><span class="hl std">,</span> <span class="hl kwc">quote</span> <span class="hl std">=</span> <span class="hl str">&quot;&quot;</span><span class="hl std">)</span>
</pre></div>
</div></div>

<p><i>df</i> is a data frame with two columns:</p>

<div class="chunk" id="unnamed-chunk-3"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl kwd">head</span><span class="hl std">(df)</span>
</pre></div>
<div class="output"><pre class="knitr r">##         log2FC gene_length_kb
## 4   0.13645050          0.019
## 10 -0.04667352          0.019
## 13 -0.07517732          0.019
## 15 -0.14915979          0.019
## 20  0.29137589          0.020
## 23 -0.09182034          0.020
</pre></div>
</div></div>

<p>it has 21334 rows:</p>

<div class="chunk" id="unnamed-chunk-4"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl kwd">dim</span><span class="hl std">(df)</span>
</pre></div>
<div class="output"><pre class="knitr r">## [1] 21334     2
</pre></div>
</div></div>

<p>To run the function and assign the results to a variable:</p>

<div class="chunk" id="unnamed-chunk-5"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl std">result</span> <span class="hl kwb">&lt;-</span> <span class="hl kwd">movingAverage.AM</span><span class="hl std">(</span><span class="hl kwc">df</span> <span class="hl std">= df ,</span> <span class="hl kwc">winSize</span> <span class="hl std">=</span> <span class="hl num">200</span> <span class="hl std">,</span> <span class="hl kwc">step</span> <span class="hl std">=</span> <span class="hl num">40</span><span class="hl std">)</span>
</pre></div>
</div></div>

<p>To extract the output table:</p>

<div class="chunk" id="unnamed-chunk-6"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl com"># result$MovingAverage</span>
<span class="hl kwd">head</span><span class="hl std">(result</span><span class="hl opt">$</span><span class="hl std">MovingAverage)</span>
</pre></div>
<div class="output"><pre class="knitr r">##   predMean      respAvg           ucl          lcl
## 1 0.028355  0.025392100  0.0602557922 -0.009471592
## 2 0.031020  0.024814286  0.0603577014 -0.010729130
## 3 0.033785  0.002299482  0.0374499663 -0.032851003
## 4 0.037580 -0.017463724  0.0159893048 -0.050916754
## 5 0.042830 -0.029813453 -0.0008006119 -0.058826294
## 6 0.049185 -0.024630348  0.0034561658 -0.052716862
</pre></div>
</div></div>

<p>To plot the output type in:</p>

<div class="chunk" id="unnamed-chunk-7"><div class="rcode"><div class="source"><pre class="knitr r"><span class="hl std">result</span><span class="hl opt">$</span><span class="hl std">plot</span>
</pre></div>
</div><div class="rimage default"><img src="figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" class="plot" /></div></div>
<p></p>
<p>GitHub: <a href="https://github.com/awnimo">awnimo</a></p>
</body>
</html>
