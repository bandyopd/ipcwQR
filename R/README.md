

<!DOCTYPE html>

<html>

<head>

<meta charset="utf-8" />
<meta name="generator" content="pandoc" />
<meta http-equiv="X-UA-Compatible" content="IE=EDGE" />

<meta name="viewport" content="width=device-width, initial-scale=1" />



</head>

<body>




<h1 class="title toc-ignore">Introduction to ipcwQR package</h1>



<div id="introduction" class="section level2">
<h2>Introduction</h2>
<p><code>ipcwQR</code> is the R package to fit the linear quantile
regressions when the data are partially interval-censored and possibly
correlated within same cluster. Let <span class="math inline">\(T\)</span> and <span class="math inline">\(X\)</span> be the event time of interest and its
related <span class="math inline">\(p\)</span>-vector of covariates,
respectively. Our main objective is to estimate the <span class="math inline">\(p\)</span>-dimensional quantile coefficient vector
<span class="math inline">\({\boldsymbol{\beta}}_0(\tau)\)</span> for
some <span class="math inline">\(\tau \in[\tau_L,\tau_R]\subset (0,
1)\)</span> in the following linear quantile regression model: <span class="math display">\[
T_i = {\bf x}_i^T {\boldsymbol{\beta}}_0(\tau) + e_i(\tau),\quad i=1,
\ldots ,n,
\]</span> where <span class="math inline">\(e_i(\tau)\)</span> is the
random error whose <span class="math inline">\(\tau\)</span>th quantile
conditional on <span class="math inline">\({\bf x}_i\)</span> equals 0.
When the data are subject to partially interval-censoring, left and
right endpoints of the censoring time, <span class="math inline">\(L\)</span> and <span class="math inline">\(R\)</span>, are observed instead of <span class="math inline">\(T\)</span> such that <span class="math inline">\(T\in(L,R)\)</span>. Note that double-censoring can
also be viewed as a special case of partly interval-censoring, i.e.,
<span class="math inline">\(T\)</span> is left-censored if <span class="math inline">\(L=0\)</span> and right-censored if <span class="math inline">\(R=\infty\)</span>.</p>
</div>
<div id="usages" class="section level2">
<h2>Usages</h2>
<p>Installation of ipcwQR package can be done by</p>
<div class="sourceCode" id="cb1"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a>devtools<span class="sc">::</span><span class="fu">install_github</span>(<span class="at">repo=</span><span class="st">&quot;YejiStat/ipcwQR&quot;</span>)</span>
<span id="cb1-2"><a href="#cb1-2" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Downloading GitHub repo YejiStat/ipcwQR@HEAD</span></span>
<span id="cb1-3"><a href="#cb1-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   </span></span>
<span id="cb1-4"><a href="#cb1-4" aria-hidden="true" tabindex="-1"></a>   checking <span class="cf">for</span> file ‘<span class="sc">/</span>private<span class="sc">/</span>var<span class="sc">/</span>folders<span class="sc">/</span>kw<span class="sc">/</span>cjk_j16d44bb_tmtfr_dpqk40000gp<span class="sc">/</span>T<span class="sc">/</span>RtmpO5WPEy<span class="sc">/</span>remotes44d73d2f736c<span class="sc">/</span>YejiStat<span class="sc">-</span>ipcwQR<span class="sc">-</span>2c2a1f2<span class="sc">/</span>DESCRIPTION’ ...</span>
<span id="cb1-5"><a href="#cb1-5" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-6"><a href="#cb1-6" aria-hidden="true" tabindex="-1"></a>✔  checking <span class="cf">for</span> file ‘<span class="sc">/</span>private<span class="sc">/</span>var<span class="sc">/</span>folders<span class="sc">/</span>kw<span class="sc">/</span>cjk_j16d44bb_tmtfr_dpqk40000gp<span class="sc">/</span>T<span class="sc">/</span>RtmpO5WPEy<span class="sc">/</span>remotes44d73d2f736c<span class="sc">/</span>YejiStat<span class="sc">-</span>ipcwQR<span class="sc">-</span>2c2a1f2<span class="sc">/</span>DESCRIPTION’</span>
<span id="cb1-7"><a href="#cb1-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-8"><a href="#cb1-8" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-9"><a href="#cb1-9" aria-hidden="true" tabindex="-1"></a>─  preparing ‘ipcwQR’<span class="sc">:</span></span>
<span id="cb1-10"><a href="#cb1-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-11"><a href="#cb1-11" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-12"><a href="#cb1-12" aria-hidden="true" tabindex="-1"></a>   checking DESCRIPTION meta<span class="sc">-</span>information ...</span>
<span id="cb1-13"><a href="#cb1-13" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-14"><a href="#cb1-14" aria-hidden="true" tabindex="-1"></a>✔  checking DESCRIPTION meta<span class="sc">-</span>information</span>
<span id="cb1-15"><a href="#cb1-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-16"><a href="#cb1-16" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-17"><a href="#cb1-17" aria-hidden="true" tabindex="-1"></a>─  checking <span class="cf">for</span> LF line<span class="sc">-</span>endings <span class="cf">in</span> source and make files and shell scripts</span>
<span id="cb1-18"><a href="#cb1-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-19"><a href="#cb1-19" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-20"><a href="#cb1-20" aria-hidden="true" tabindex="-1"></a>─  checking <span class="cf">for</span> empty or unneeded directories</span>
<span id="cb1-21"><a href="#cb1-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-22"><a href="#cb1-22" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-23"><a href="#cb1-23" aria-hidden="true" tabindex="-1"></a>─  building ‘ipcwQR_0.<span class="dv">0</span>.<span class="dv">0</span>.<span class="fl">9000.</span>tar.gz’</span>
<span id="cb1-24"><a href="#cb1-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-25"><a href="#cb1-25" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb1-26"><a href="#cb1-26" aria-hidden="true" tabindex="-1"></a>   </span>
<span id="cb1-27"><a href="#cb1-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb1-28"><a href="#cb1-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Warning in i.p(...): installation of package &#39;/var/folders/kw/</span></span>
<span id="cb1-29"><a href="#cb1-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; cjk_j16d44bb_tmtfr_dpqk40000gp/T//RtmpO5WPEy/file44d77642d6d5/</span></span>
<span id="cb1-30"><a href="#cb1-30" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ipcwQR_0.0.0.9000.tar.gz&#39; had non-zero exit status</span></span></code></pre></div>
<p>or</p>
<div class="sourceCode" id="cb2"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a>base<span class="sc">::</span><span class="fu">require</span>(<span class="st">&quot;ipcwQR&quot;</span>)</span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: ipcwQR</span></span></code></pre></div>
<p>picrq() function has the following arguments:</p>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="co"># picrq(L=V,R=U,delta=delta,x=x,tau=tau)</span></span></code></pre></div>
<p>see the detailed description from help(picrq()).</p>
<p>We first simulate univariate partly interval-censored data with
normal random error, which is a similar simulation setting of Kim et
al. (2022+).</p>
<div class="sourceCode" id="cb4"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb4-1"><a href="#cb4-1" aria-hidden="true" tabindex="-1"></a>PICdata<span class="ot">=</span><span class="cf">function</span>(n,tau,case){</span>
<span id="cb4-2"><a href="#cb4-2" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb4-3"><a href="#cb4-3" aria-hidden="true" tabindex="-1"></a>  <span class="fu">library</span>(tidyverse)</span>
<span id="cb4-4"><a href="#cb4-4" aria-hidden="true" tabindex="-1"></a>  <span class="fu">library</span>(survival)</span>
<span id="cb4-5"><a href="#cb4-5" aria-hidden="true" tabindex="-1"></a>  <span class="cf">if</span> (case<span class="sc">==</span><span class="st">&quot;norm85&quot;</span>){err <span class="ot">=</span> <span class="fu">rnorm</span>(n); q<span class="ot">=</span><span class="fu">qnorm</span>(tau); p0 <span class="ot">=</span> <span class="fl">0.90</span>}</span>
<span id="cb4-6"><a href="#cb4-6" aria-hidden="true" tabindex="-1"></a>  <span class="cf">else</span> <span class="cf">if</span> (case<span class="sc">==</span><span class="st">&quot;norm75&quot;</span>){err <span class="ot">=</span> <span class="fu">rnorm</span>(n); q<span class="ot">=</span><span class="fu">qnorm</span>(tau); p0 <span class="ot">=</span> <span class="fl">0.85</span>}</span>
<span id="cb4-7"><a href="#cb4-7" aria-hidden="true" tabindex="-1"></a>  <span class="cf">else</span> <span class="cf">if</span> (case<span class="sc">==</span><span class="st">&quot;gum85&quot;</span>){err<span class="ot">=</span><span class="fu">revd</span>(n); q<span class="ot">=</span><span class="fu">qevd</span>(tau); p0 <span class="ot">=</span> <span class="fl">0.91</span>}</span>
<span id="cb4-8"><a href="#cb4-8" aria-hidden="true" tabindex="-1"></a>  <span class="cf">else</span> <span class="cf">if</span> (case<span class="sc">==</span><span class="st">&quot;gum75&quot;</span>){err<span class="ot">=</span><span class="fu">revd</span>(n); q<span class="ot">=</span><span class="fu">qevd</span>(tau);  p0 <span class="ot">=</span> <span class="fl">0.84</span>}</span>
<span id="cb4-9"><a href="#cb4-9" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb4-10"><a href="#cb4-10" aria-hidden="true" tabindex="-1"></a>  x1<span class="ot">=</span><span class="fu">runif</span>(n,<span class="sc">-</span><span class="fl">1.2</span>,<span class="fl">1.7</span>); x2<span class="ot">=</span><span class="fu">rbinom</span>(n,<span class="dv">1</span>,<span class="fl">0.6</span>)</span>
<span id="cb4-11"><a href="#cb4-11" aria-hidden="true" tabindex="-1"></a>  Time <span class="ot">=</span> <span class="fl">1.7</span><span class="sc">+</span>x1<span class="sc">+</span>x2<span class="sc">+</span>err<span class="sc">*</span>(<span class="dv">1</span><span class="fl">-0.1</span><span class="sc">*</span>x2)</span>
<span id="cb4-12"><a href="#cb4-12" aria-hidden="true" tabindex="-1"></a>  int <span class="ot">=</span> <span class="fu">outer</span>(<span class="dv">1</span><span class="sc">:</span>n, <span class="fu">c</span>(<span class="dv">0</span>,<span class="dv">0</span>), <span class="st">&quot;*&quot;</span>)</span>
<span id="cb4-13"><a href="#cb4-13" aria-hidden="true" tabindex="-1"></a>  delta <span class="ot">=</span> <span class="dv">0</span></span>
<span id="cb4-14"><a href="#cb4-14" aria-hidden="true" tabindex="-1"></a>  prob <span class="ot">=</span> p0 <span class="sc">-</span> <span class="fl">0.05</span><span class="sc">*</span>(x1 <span class="sc">&gt;</span> <span class="dv">0</span>) <span class="co">#950 963 gum75%0.5</span></span>
<span id="cb4-15"><a href="#cb4-15" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb4-16"><a href="#cb4-16" aria-hidden="true" tabindex="-1"></a>  <span class="cf">for</span> (i <span class="cf">in</span> <span class="dv">1</span><span class="sc">:</span>n) {</span>
<span id="cb4-17"><a href="#cb4-17" aria-hidden="true" tabindex="-1"></a>    ind <span class="ot">=</span> <span class="fu">rbinom</span>(<span class="dv">1</span>, <span class="dv">1</span>, prob[i])</span>
<span id="cb4-18"><a href="#cb4-18" aria-hidden="true" tabindex="-1"></a>    <span class="cf">if</span> (ind <span class="sc">==</span> <span class="dv">1</span>) {</span>
<span id="cb4-19"><a href="#cb4-19" aria-hidden="true" tabindex="-1"></a>      int[i,] <span class="ot">=</span> <span class="fu">rep</span>(Time[i], <span class="dv">2</span>)</span>
<span id="cb4-20"><a href="#cb4-20" aria-hidden="true" tabindex="-1"></a>      delta[i] <span class="ot">=</span> <span class="dv">1</span></span>
<span id="cb4-21"><a href="#cb4-21" aria-hidden="true" tabindex="-1"></a>    } <span class="cf">else</span> {</span>
<span id="cb4-22"><a href="#cb4-22" aria-hidden="true" tabindex="-1"></a>      utmp <span class="ot">=</span> U <span class="ot">=</span> <span class="fl">1e-8</span></span>
<span id="cb4-23"><a href="#cb4-23" aria-hidden="true" tabindex="-1"></a>      j<span class="ot">=</span><span class="dv">1</span></span>
<span id="cb4-24"><a href="#cb4-24" aria-hidden="true" tabindex="-1"></a>      <span class="cf">while</span> (utmp <span class="sc">&lt;=</span> <span class="dv">300</span>) {</span>
<span id="cb4-25"><a href="#cb4-25" aria-hidden="true" tabindex="-1"></a>        <span class="co"># U[j+1] = U[j] + runif(1,0.1,1)</span></span>
<span id="cb4-26"><a href="#cb4-26" aria-hidden="true" tabindex="-1"></a>        U[j<span class="sc">+</span><span class="dv">1</span>] <span class="ot">=</span> U[j] <span class="sc">+</span> <span class="fu">rexp</span>(<span class="dv">1</span>,<span class="dv">1</span>)</span>
<span id="cb4-27"><a href="#cb4-27" aria-hidden="true" tabindex="-1"></a>        utmp <span class="ot">=</span> U[j<span class="sc">+</span><span class="dv">1</span>]</span>
<span id="cb4-28"><a href="#cb4-28" aria-hidden="true" tabindex="-1"></a>        j <span class="ot">=</span> j<span class="sc">+</span><span class="dv">1</span></span>
<span id="cb4-29"><a href="#cb4-29" aria-hidden="true" tabindex="-1"></a>      }</span>
<span id="cb4-30"><a href="#cb4-30" aria-hidden="true" tabindex="-1"></a>      L <span class="ot">=</span> U[<span class="dv">1</span><span class="sc">:</span>(j<span class="dv">-2</span>)]</span>
<span id="cb4-31"><a href="#cb4-31" aria-hidden="true" tabindex="-1"></a>      R <span class="ot">=</span> U[<span class="dv">2</span><span class="sc">:</span>(j<span class="dv">-1</span>)]</span>
<span id="cb4-32"><a href="#cb4-32" aria-hidden="true" tabindex="-1"></a>      </span>
<span id="cb4-33"><a href="#cb4-33" aria-hidden="true" tabindex="-1"></a>      <span class="cf">if</span> (Time[i] <span class="sc">&lt;</span> <span class="fu">min</span>(L)) {</span>
<span id="cb4-34"><a href="#cb4-34" aria-hidden="true" tabindex="-1"></a>        int[i,] <span class="ot">=</span>  <span class="fu">c</span>(<span class="fl">1e-8</span>, <span class="fu">min</span>(L))</span>
<span id="cb4-35"><a href="#cb4-35" aria-hidden="true" tabindex="-1"></a>        delta[i] <span class="ot">=</span> <span class="dv">0</span></span>
<span id="cb4-36"><a href="#cb4-36" aria-hidden="true" tabindex="-1"></a>      } <span class="cf">else</span> <span class="cf">if</span> (Time[i] <span class="sc">&gt;</span> <span class="fu">max</span>(R)) {</span>
<span id="cb4-37"><a href="#cb4-37" aria-hidden="true" tabindex="-1"></a>        int[i,] <span class="ot">=</span> <span class="fu">c</span>(<span class="fu">max</span>(R), <span class="fl">1e8</span>)</span>
<span id="cb4-38"><a href="#cb4-38" aria-hidden="true" tabindex="-1"></a>        delta[i] <span class="ot">=</span> <span class="dv">0</span></span>
<span id="cb4-39"><a href="#cb4-39" aria-hidden="true" tabindex="-1"></a>      } <span class="cf">else</span> {</span>
<span id="cb4-40"><a href="#cb4-40" aria-hidden="true" tabindex="-1"></a>        idd <span class="ot">=</span> (Time[i] <span class="sc">&gt;</span> L)<span class="sc">&amp;</span>(Time[i] <span class="sc">&lt;</span> R)</span>
<span id="cb4-41"><a href="#cb4-41" aria-hidden="true" tabindex="-1"></a>        <span class="cf">if</span> (<span class="fu">sum</span>(idd) <span class="sc">==</span> <span class="dv">1</span>) {</span>
<span id="cb4-42"><a href="#cb4-42" aria-hidden="true" tabindex="-1"></a>          int[i,] <span class="ot">=</span> <span class="fu">t</span>(<span class="fu">cbind</span>(L,R)[idd,])</span>
<span id="cb4-43"><a href="#cb4-43" aria-hidden="true" tabindex="-1"></a>          delta[i] <span class="ot">=</span> <span class="dv">0</span></span>
<span id="cb4-44"><a href="#cb4-44" aria-hidden="true" tabindex="-1"></a>        }</span>
<span id="cb4-45"><a href="#cb4-45" aria-hidden="true" tabindex="-1"></a>        <span class="cf">if</span> (<span class="fu">sum</span>(idd) <span class="sc">==</span> <span class="dv">0</span>) {</span>
<span id="cb4-46"><a href="#cb4-46" aria-hidden="true" tabindex="-1"></a>          int[i,] <span class="ot">=</span> <span class="fu">rep</span>(Time[i], <span class="dv">2</span>)</span>
<span id="cb4-47"><a href="#cb4-47" aria-hidden="true" tabindex="-1"></a>          delta[i] <span class="ot">=</span> <span class="dv">1</span></span>
<span id="cb4-48"><a href="#cb4-48" aria-hidden="true" tabindex="-1"></a>        }</span>
<span id="cb4-49"><a href="#cb4-49" aria-hidden="true" tabindex="-1"></a>      }</span>
<span id="cb4-50"><a href="#cb4-50" aria-hidden="true" tabindex="-1"></a>    }</span>
<span id="cb4-51"><a href="#cb4-51" aria-hidden="true" tabindex="-1"></a>  }</span>
<span id="cb4-52"><a href="#cb4-52" aria-hidden="true" tabindex="-1"></a>  dt <span class="ot">=</span> <span class="fu">data.frame</span>(int, delta, x1, x2)</span>
<span id="cb4-53"><a href="#cb4-53" aria-hidden="true" tabindex="-1"></a>  <span class="fu">colnames</span>(dt) <span class="ot">=</span> <span class="fu">c</span>(<span class="st">&quot;L&quot;</span>, <span class="st">&quot;R&quot;</span>, <span class="st">&quot;delta&quot;</span>, <span class="st">&quot;x1&quot;</span>, <span class="st">&quot;x2&quot;</span>)</span>
<span id="cb4-54"><a href="#cb4-54" aria-hidden="true" tabindex="-1"></a>  dt</span>
<span id="cb4-55"><a href="#cb4-55" aria-hidden="true" tabindex="-1"></a>}</span>
<span id="cb4-56"><a href="#cb4-56" aria-hidden="true" tabindex="-1"></a>n<span class="ot">=</span><span class="dv">200</span>; case<span class="ot">=</span><span class="st">&quot;norm85&quot;</span>; tau<span class="ot">=</span><span class="fl">0.3</span></span>
<span id="cb4-57"><a href="#cb4-57" aria-hidden="true" tabindex="-1"></a>d<span class="ot">=</span><span class="fu">PICdata</span>(<span class="at">n=</span>n,<span class="at">case =</span> case,<span class="at">tau=</span>tau)</span>
<span id="cb4-58"><a href="#cb4-58" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──</span></span>
<span id="cb4-59"><a href="#cb4-59" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 </span></span>
<span id="cb4-60"><a href="#cb4-60" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ tibble  3.1.8      ✔ dplyr   1.0.10</span></span>
<span id="cb4-61"><a href="#cb4-61" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ tidyr   1.2.0      ✔ stringr 1.4.0 </span></span>
<span id="cb4-62"><a href="#cb4-62" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ readr   2.1.2      ✔ forcats 0.5.2 </span></span>
<span id="cb4-63"><a href="#cb4-63" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──</span></span>
<span id="cb4-64"><a href="#cb4-64" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✖ dplyr::filter() masks stats::filter()</span></span>
<span id="cb4-65"><a href="#cb4-65" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✖ dplyr::lag()    masks stats::lag()</span></span>
<span id="cb4-66"><a href="#cb4-66" aria-hidden="true" tabindex="-1"></a>x<span class="ot">=</span><span class="fu">with</span>(d, <span class="fu">cbind</span>(x1,x2))</span>
<span id="cb4-67"><a href="#cb4-67" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(<span class="at">L=</span>d<span class="sc">$</span>L,<span class="at">R=</span>d<span class="sc">$</span>R,<span class="at">delta=</span>d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau)</span>
<span id="cb4-68"><a href="#cb4-68" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb4-69"><a href="#cb4-69" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.444219 0.172798  0e+00 1.105534 1.782904</span></span>
<span id="cb4-70"><a href="#cb4-70" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     1.039211 0.124371  0e+00 0.795444 1.282979</span></span>
<span id="cb4-71"><a href="#cb4-71" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     0.925415 0.225434  2e-05 0.483565 1.367264</span></span>
<span id="cb4-72"><a href="#cb4-72" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(<span class="at">L=</span>d<span class="sc">$</span>L,<span class="at">R=</span>d<span class="sc">$</span>R,<span class="at">delta=</span>d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">estimation =</span> <span class="st">&quot;dr&quot;</span>)</span>
<span id="cb4-73"><a href="#cb4-73" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb4-74"><a href="#cb4-74" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.444172 0.172807  0e+00 1.105470 1.782873</span></span>
<span id="cb4-75"><a href="#cb4-75" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     1.039211 0.124377  0e+00 0.795433 1.282989</span></span>
<span id="cb4-76"><a href="#cb4-76" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     0.925415 0.225443  2e-05 0.483545 1.367284</span></span></code></pre></div>
<p>We posit two estimating methods, ipcw method and doubly robust ipcw
method, which can be conducted by specifying estimation = NULL and
estimation = “dr”, respectively.</p>
<div class="sourceCode" id="cb5"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a>DCdata <span class="ot">=</span> <span class="cf">function</span>(n, <span class="at">case=</span>case){</span>
<span id="cb5-2"><a href="#cb5-2" aria-hidden="true" tabindex="-1"></a>  <span class="cf">if</span>(case<span class="sc">==</span><span class="st">&quot;norm85&quot;</span>){err<span class="ot">=</span><span class="fu">rnorm</span>(n); l<span class="ot">=</span><span class="fl">1.9</span>; r<span class="ot">=</span><span class="fl">8.1</span></span>
<span id="cb5-3"><a href="#cb5-3" aria-hidden="true" tabindex="-1"></a>  }<span class="cf">else</span> <span class="cf">if</span>(case<span class="sc">==</span><span class="st">&quot;norm75&quot;</span>){err<span class="ot">=</span><span class="fu">rnorm</span>(n); l<span class="ot">=</span><span class="fl">2.7</span>; r<span class="ot">=</span><span class="fl">6.2</span></span>
<span id="cb5-4"><a href="#cb5-4" aria-hidden="true" tabindex="-1"></a>  }<span class="cf">else</span> <span class="cf">if</span>(case<span class="sc">==</span><span class="st">&quot;gum85&quot;</span>){err<span class="ot">=</span><span class="fu">revd</span>(n);  l<span class="ot">=</span><span class="fl">2.3</span>; r<span class="ot">=</span><span class="fl">11.5</span></span>
<span id="cb5-5"><a href="#cb5-5" aria-hidden="true" tabindex="-1"></a>  }<span class="cf">else</span> <span class="cf">if</span>(case<span class="sc">==</span><span class="st">&quot;gum75&quot;</span>){err<span class="ot">=</span><span class="fu">revd</span>(n); l<span class="ot">=</span><span class="fl">3.3</span>; r<span class="ot">=</span><span class="fl">7.9</span></span>
<span id="cb5-6"><a href="#cb5-6" aria-hidden="true" tabindex="-1"></a>  }</span>
<span id="cb5-7"><a href="#cb5-7" aria-hidden="true" tabindex="-1"></a>  x1<span class="ot">=</span><span class="fu">runif</span>(n,<span class="sc">-</span><span class="fl">1.2</span>,<span class="fl">1.7</span>); x2<span class="ot">=</span><span class="fu">rbinom</span>(n,<span class="dv">1</span>,<span class="fl">0.6</span>)</span>
<span id="cb5-8"><a href="#cb5-8" aria-hidden="true" tabindex="-1"></a>  T <span class="ot">=</span> <span class="fl">1.7</span><span class="sc">+</span>x1<span class="sc">+</span>x2<span class="sc">+</span>err<span class="sc">*</span>(<span class="dv">1</span><span class="fl">-0.1</span><span class="sc">*</span>x2)</span>
<span id="cb5-9"><a href="#cb5-9" aria-hidden="true" tabindex="-1"></a>  L<span class="ot">=</span><span class="fu">runif</span>(n,<span class="sc">-</span><span class="fl">2.8</span>,l); R<span class="ot">=</span>L<span class="sc">+</span><span class="fu">runif</span>(n,<span class="fl">4.2</span>,r)</span>
<span id="cb5-10"><a href="#cb5-10" aria-hidden="true" tabindex="-1"></a>  Y<span class="ot">=</span><span class="fu">pmin</span>(R,<span class="fu">pmax</span>(T,L))</span>
<span id="cb5-11"><a href="#cb5-11" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb5-12"><a href="#cb5-12" aria-hidden="true" tabindex="-1"></a>  delta<span class="ot">=</span><span class="fu">case_when</span>(</span>
<span id="cb5-13"><a href="#cb5-13" aria-hidden="true" tabindex="-1"></a>    T<span class="sc">&lt;</span>L <span class="sc">~</span> <span class="dv">3</span>,</span>
<span id="cb5-14"><a href="#cb5-14" aria-hidden="true" tabindex="-1"></a>    T<span class="sc">&gt;</span>R <span class="sc">~</span> <span class="dv">2</span>,</span>
<span id="cb5-15"><a href="#cb5-15" aria-hidden="true" tabindex="-1"></a>    <span class="cn">TRUE</span> <span class="sc">~</span> <span class="dv">1</span> <span class="co">#observed</span></span>
<span id="cb5-16"><a href="#cb5-16" aria-hidden="true" tabindex="-1"></a>  )</span>
<span id="cb5-17"><a href="#cb5-17" aria-hidden="true" tabindex="-1"></a>  <span class="fu">table</span>(delta)<span class="sc">/</span><span class="fu">length</span>(Y)</span>
<span id="cb5-18"><a href="#cb5-18" aria-hidden="true" tabindex="-1"></a>  d<span class="ot">=</span><span class="fu">data.frame</span>(<span class="at">Y=</span>Y,<span class="at">L=</span>L,<span class="at">R=</span>R,<span class="at">T=</span>T,<span class="at">delta=</span>delta,<span class="at">x1=</span>x1,<span class="at">x2=</span>x2)</span>
<span id="cb5-19"><a href="#cb5-19" aria-hidden="true" tabindex="-1"></a>  d[<span class="fu">order</span>(Y),]</span>
<span id="cb5-20"><a href="#cb5-20" aria-hidden="true" tabindex="-1"></a>}</span>
<span id="cb5-21"><a href="#cb5-21" aria-hidden="true" tabindex="-1"></a>n<span class="ot">=</span><span class="dv">200</span>; case<span class="ot">=</span><span class="st">&quot;norm85&quot;</span>; tau<span class="ot">=</span><span class="fl">0.3</span></span>
<span id="cb5-22"><a href="#cb5-22" aria-hidden="true" tabindex="-1"></a>d<span class="ot">=</span><span class="fu">DCdata</span>(<span class="at">n=</span>n,<span class="at">case =</span> case)</span>
<span id="cb5-23"><a href="#cb5-23" aria-hidden="true" tabindex="-1"></a>x<span class="ot">=</span><span class="fu">with</span>(d, <span class="fu">cbind</span>(x1,x2))</span>
<span id="cb5-24"><a href="#cb5-24" aria-hidden="true" tabindex="-1"></a>L<span class="ot">=</span>d<span class="sc">$</span>L; R<span class="ot">=</span>d<span class="sc">$</span>R; T<span class="ot">=</span>d<span class="sc">$</span>T; delta<span class="ot">=</span>d<span class="sc">$</span>delta; x<span class="ot">=</span><span class="fu">with</span>(d, <span class="fu">cbind</span>(x1,x2))</span>
<span id="cb5-25"><a href="#cb5-25" aria-hidden="true" tabindex="-1"></a><span class="fu">dcrq</span>(L,R,T,delta,x,tau)</span>
<span id="cb5-26"><a href="#cb5-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb5-27"><a href="#cb5-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.232325 0.241320  0e+00 0.759338 1.705312</span></span>
<span id="cb5-28"><a href="#cb5-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     0.774301 0.145372  0e+00 0.489373 1.059230</span></span>
<span id="cb5-29"><a href="#cb5-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     1.240393 0.276294  4e-06 0.698856 1.781929</span></span>
<span id="cb5-30"><a href="#cb5-30" aria-hidden="true" tabindex="-1"></a><span class="fu">dcrq</span>(L,R,T,delta,x,tau,<span class="at">estimation =</span> <span class="st">&quot;dr&quot;</span>)</span>
<span id="cb5-31"><a href="#cb5-31" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb5-32"><a href="#cb5-32" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.227054 0.242834  0e+00 0.751099 1.703008</span></span>
<span id="cb5-33"><a href="#cb5-33" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     0.774745 0.146693  0e+00 0.487228 1.062263</span></span>
<span id="cb5-34"><a href="#cb5-34" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     1.239997 0.278410  4e-06 0.694314 1.785680</span></span></code></pre></div>
<p>Next, we give illustrative example of the method for the multivariate
clustered data, by using phase 3 metastatic colorectal cancer clinical
trial. This dataset is available data(mCRC) in PICBayes package (Pan,
2021).</p>
<div class="sourceCode" id="cb6"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(PICBayes)</span>
<span id="cb6-2"><a href="#cb6-2" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: coda</span></span>
<span id="cb6-3"><a href="#cb6-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: MCMCpack</span></span>
<span id="cb6-4"><a href="#cb6-4" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: MASS</span></span>
<span id="cb6-5"><a href="#cb6-5" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb6-6"><a href="#cb6-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Attaching package: &#39;MASS&#39;</span></span>
<span id="cb6-7"><a href="#cb6-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; The following object is masked from &#39;package:dplyr&#39;:</span></span>
<span id="cb6-8"><a href="#cb6-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb6-9"><a href="#cb6-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     select</span></span>
<span id="cb6-10"><a href="#cb6-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb6-11"><a href="#cb6-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Markov Chain Monte Carlo Package (MCMCpack)</span></span>
<span id="cb6-12"><a href="#cb6-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Copyright (C) 2003-2022 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park</span></span>
<span id="cb6-13"><a href="#cb6-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb6-14"><a href="#cb6-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Support provided by the U.S. National Science Foundation</span></span>
<span id="cb6-15"><a href="#cb6-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## (Grants SES-0350646 and SES-0350613)</span></span>
<span id="cb6-16"><a href="#cb6-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb6-17"><a href="#cb6-17" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-18"><a href="#cb6-18" aria-hidden="true" tabindex="-1"></a><span class="fu">data</span>(<span class="st">&quot;mCRC&quot;</span>)</span>
<span id="cb6-19"><a href="#cb6-19" aria-hidden="true" tabindex="-1"></a>d <span class="ot">=</span> <span class="fu">with</span>(<span class="fu">data.frame</span>(mCRC), <span class="fu">data.frame</span>(<span class="at">L =</span> <span class="fu">as.numeric</span>(L),</span>
<span id="cb6-20"><a href="#cb6-20" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">R =</span> <span class="fu">as.numeric</span>(R),</span>
<span id="cb6-21"><a href="#cb6-21" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">U =</span> <span class="fu">ifelse</span>(y<span class="sc">==</span><span class="dv">0</span>,R,L),</span>
<span id="cb6-22"><a href="#cb6-22" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">V =</span> <span class="fu">ifelse</span>(y<span class="sc">==</span><span class="dv">2</span>,L,R),</span>
<span id="cb6-23"><a href="#cb6-23" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Cluster weighted data</span></span>
<span id="cb6-24"><a href="#cb6-24" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">id=</span>(<span class="fu">rep</span>(<span class="fu">c</span>(<span class="fu">table</span>(SITE)),<span class="fu">c</span>(<span class="fu">table</span>(SITE)))),</span>
<span id="cb6-25"><a href="#cb6-25" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.</span></span>
<span id="cb6-26"><a href="#cb6-26" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">x1=</span> <span class="fu">case_when</span>(TRT_C <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">0</span>, <span class="co">#Pan et al data</span></span>
<span id="cb6-27"><a href="#cb6-27" aria-hidden="true" tabindex="-1"></a>                                                    TRT_C <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">1</span>),</span>
<span id="cb6-28"><a href="#cb6-28" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.</span></span>
<span id="cb6-29"><a href="#cb6-29" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">x2=</span> <span class="fu">case_when</span>(KRAS_C <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">1</span>,</span>
<span id="cb6-30"><a href="#cb6-30" aria-hidden="true" tabindex="-1"></a>                                                    KRAS_C <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">0</span>),</span>
<span id="cb6-31"><a href="#cb6-31" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">site =</span> <span class="fu">as.numeric</span>(SITE),</span>
<span id="cb6-32"><a href="#cb6-32" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">y =</span> <span class="fu">as.numeric</span>(y),</span>
<span id="cb6-33"><a href="#cb6-33" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">delta =</span> <span class="fu">case_when</span>(IC <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">1</span>,</span>
<span id="cb6-34"><a href="#cb6-34" aria-hidden="true" tabindex="-1"></a>                                                        IC <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">0</span>)</span>
<span id="cb6-35"><a href="#cb6-35" aria-hidden="true" tabindex="-1"></a>));</span>
<span id="cb6-36"><a href="#cb6-36" aria-hidden="true" tabindex="-1"></a>L<span class="ot">=</span>d<span class="sc">$</span>U;R<span class="ot">=</span>d<span class="sc">$</span>V; delta<span class="ot">=</span>d<span class="sc">$</span>delta</span>
<span id="cb6-37"><a href="#cb6-37" aria-hidden="true" tabindex="-1"></a>L<span class="ot">=</span>(<span class="fu">log</span>(d<span class="sc">$</span>U));R<span class="ot">=</span><span class="fu">log</span>(d<span class="sc">$</span>V); delta<span class="ot">=</span>d<span class="sc">$</span>delta</span>
<span id="cb6-38"><a href="#cb6-38" aria-hidden="true" tabindex="-1"></a>x <span class="ot">=</span> <span class="fu">cbind</span>(d<span class="sc">$</span>x1,d<span class="sc">$</span>x2); id<span class="ot">=</span>d<span class="sc">$</span>id;  tau<span class="ot">=</span><span class="fl">0.3</span>;</span>
<span id="cb6-39"><a href="#cb6-39" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-40"><a href="#cb6-40" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(L,R,delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau)</span>
<span id="cb6-41"><a href="#cb6-41" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb6-42"><a href="#cb6-42" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     4.288892 0.535586 0.000000  3.239142 5.338641</span></span>
<span id="cb6-43"><a href="#cb6-43" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.144972 0.552113 0.396438 -0.937169 1.227113</span></span>
<span id="cb6-44"><a href="#cb6-44" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.425358 0.599873 0.239137 -0.750392 1.601109</span></span>
<span id="cb6-45"><a href="#cb6-45" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(L,R,delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau, <span class="at">estimation =</span> <span class="st">&quot;dr&quot;</span>)</span>
<span id="cb6-46"><a href="#cb6-46" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb6-47"><a href="#cb6-47" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     4.288876 0.535590 0.000000  3.239119 5.338632</span></span>
<span id="cb6-48"><a href="#cb6-48" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.144972 0.552116 0.396439 -0.937174 1.227118</span></span>
<span id="cb6-49"><a href="#cb6-49" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.425369 0.599876 0.239133 -0.750388 1.601127</span></span>
<span id="cb6-50"><a href="#cb6-50" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(L,R,delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">id=</span>id,<span class="at">h=</span><span class="fl">0.9</span>,<span class="at">k=</span><span class="dv">2</span>)</span>
<span id="cb6-51"><a href="#cb6-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb6-52"><a href="#cb6-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     4.334544 0.012620      0 4.309808 4.359280</span></span>
<span id="cb6-53"><a href="#cb6-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.093215 0.013290      0 0.067166 0.119263</span></span>
<span id="cb6-54"><a href="#cb6-54" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.417785 0.014393      0 0.389575 0.445994</span></span>
<span id="cb6-55"><a href="#cb6-55" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(L,R,delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">id=</span><span class="cn">NULL</span>,<span class="at">h=</span><span class="fl">0.9</span>)</span>
<span id="cb6-56"><a href="#cb6-56" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb6-57"><a href="#cb6-57" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     4.288892 0.535586 0.000000  3.239142 5.338641</span></span>
<span id="cb6-58"><a href="#cb6-58" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.144972 0.552113 0.396438 -0.937169 1.227113</span></span>
<span id="cb6-59"><a href="#cb6-59" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.425358 0.599873 0.239137 -0.750392 1.601109</span></span></code></pre></div>
</div>
<div id="references" class="section level2">
<h2>References</h2>
<ul>
<li><p>Chiou, S. H., Kang, S. and Yan, J. (2015). Rank-based estimating
equations with general weight for accelerated failure time models: an
induced smoothing approach. Statistics in Medicine 34(9):
1495–-1510.</p></li>
<li><p>Pan, C. (2021). PICBayes: Bayesian Models for Partly
Interval-Censored Data. R package. <a href="https://CRAN.R-project.org/package=PICBayes" class="uri">https://CRAN.R-project.org/package=PICBayes</a>.</p></li>
<li><p>Kim, Y., Choi, T., Park, S., Choi, S. and Bandyopadhyay, D.
(2022+). Inverse weighted quantile regression with partially
interval-censored data. <em>Submitted to SMMR</em>.</p></li>
</ul>
</div>
