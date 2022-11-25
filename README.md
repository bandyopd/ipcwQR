# ipcwQR

<!DOCTYPE html>

<html>

<head>

<meta charset="utf-8" />
<meta name="generator" content="pandoc" />
<meta http-equiv="X-UA-Compatible" content="IE=EDGE" />

<meta name="viewport" content="width=device-width, initial-scale=1" />



<title>Introduction to ipcwQR package</title>

<script>// Pandoc 2.9 adds attributes on both header and div. We remove the former (to
// be compatible with the behavior of Pandoc < 2.8).
document.addEventListener('DOMContentLoaded', function(e) {
  var hs = document.querySelectorAll("div.section[class*='level'] > :first-child");
  var i, h, a;
  for (i = 0; i < hs.length; i++) {
    h = hs[i];
    if (!/^h[1-6]$/i.test(h.tagName)) continue;  // it should be a header h1-h6
    a = h.attributes;
    while (a.length > 0) h.removeAttribute(a[0].name);
  }
});
</script>

<style type="text/css">
  code{white-space: pre-wrap;}
  span.smallcaps{font-variant: small-caps;}
  span.underline{text-decoration: underline;}
  div.column{display: inline-block; vertical-align: top; width: 50%;}
  div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
  ul.task-list{list-style: none;}
    </style>



<style type="text/css">
  code {
    white-space: pre;
  }
  .sourceCode {
    overflow: visible;
  }
</style>
<style type="text/css" data-origin="pandoc">
pre > code.sourceCode { white-space: pre; position: relative; }
pre > code.sourceCode > span { display: inline-block; line-height: 1.25; }
pre > code.sourceCode > span:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode > span { color: inherit; text-decoration: inherit; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
pre > code.sourceCode { white-space: pre-wrap; }
pre > code.sourceCode > span { text-indent: -5em; padding-left: 5em; }
}
pre.numberSource code
  { counter-reset: source-line 0; }
pre.numberSource code > span
  { position: relative; left: -4em; counter-increment: source-line; }
pre.numberSource code > span > a:first-child::before
  { content: counter(source-line);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
    color: #aaaaaa;
  }
pre.numberSource { margin-left: 3em; border-left: 1px solid #aaaaaa;  padding-left: 4px; }
div.sourceCode
  {   }
@media screen {
pre > code.sourceCode > span > a:first-child::before { text-decoration: underline; }
}
code span.al { color: #ff0000; font-weight: bold; } /* Alert */
code span.an { color: #60a0b0; font-weight: bold; font-style: italic; } /* Annotation */
code span.at { color: #7d9029; } /* Attribute */
code span.bn { color: #40a070; } /* BaseN */
code span.bu { } /* BuiltIn */
code span.cf { color: #007020; font-weight: bold; } /* ControlFlow */
code span.ch { color: #4070a0; } /* Char */
code span.cn { color: #880000; } /* Constant */
code span.co { color: #60a0b0; font-style: italic; } /* Comment */
code span.cv { color: #60a0b0; font-weight: bold; font-style: italic; } /* CommentVar */
code span.do { color: #ba2121; font-style: italic; } /* Documentation */
code span.dt { color: #902000; } /* DataType */
code span.dv { color: #40a070; } /* DecVal */
code span.er { color: #ff0000; font-weight: bold; } /* Error */
code span.ex { } /* Extension */
code span.fl { color: #40a070; } /* Float */
code span.fu { color: #06287e; } /* Function */
code span.im { } /* Import */
code span.in { color: #60a0b0; font-weight: bold; font-style: italic; } /* Information */
code span.kw { color: #007020; font-weight: bold; } /* Keyword */
code span.op { color: #666666; } /* Operator */
code span.ot { color: #007020; } /* Other */
code span.pp { color: #bc7a00; } /* Preprocessor */
code span.sc { color: #4070a0; } /* SpecialChar */
code span.ss { color: #bb6688; } /* SpecialString */
code span.st { color: #4070a0; } /* String */
code span.va { color: #19177c; } /* Variable */
code span.vs { color: #4070a0; } /* VerbatimString */
code span.wa { color: #60a0b0; font-weight: bold; font-style: italic; } /* Warning */

</style>
<script>
// apply pandoc div.sourceCode style to pre.sourceCode instead
(function() {
  var sheets = document.styleSheets;
  for (var i = 0; i < sheets.length; i++) {
    if (sheets[i].ownerNode.dataset["origin"] !== "pandoc") continue;
    try { var rules = sheets[i].cssRules; } catch (e) { continue; }
    var j = 0;
    while (j < rules.length) {
      var rule = rules[j];
      // check if there is a div.sourceCode rule
      if (rule.type !== rule.STYLE_RULE || rule.selectorText !== "div.sourceCode") {
        j++;
        continue;
      }
      var style = rule.style.cssText;
      // check if color or background-color is set
      if (rule.style.color === '' && rule.style.backgroundColor === '') {
        j++;
        continue;
      }
      // replace div.sourceCode by a pre.sourceCode rule
      sheets[i].deleteRule(j);
      sheets[i].insertRule('pre.sourceCode{' + style + '}', j);
    }
  }
})();
</script>




<style type="text/css">body {
background-color: #fff;
margin: 1em auto;
max-width: 700px;
overflow: visible;
padding-left: 2em;
padding-right: 2em;
font-family: "Open Sans", "Helvetica Neue", Helvetica, Arial, sans-serif;
font-size: 14px;
line-height: 1.35;
}
#TOC {
clear: both;
margin: 0 0 10px 10px;
padding: 4px;
width: 400px;
border: 1px solid #CCCCCC;
border-radius: 5px;
background-color: #f6f6f6;
font-size: 13px;
line-height: 1.3;
}
#TOC .toctitle {
font-weight: bold;
font-size: 15px;
margin-left: 5px;
}
#TOC ul {
padding-left: 40px;
margin-left: -1.5em;
margin-top: 5px;
margin-bottom: 5px;
}
#TOC ul ul {
margin-left: -2em;
}
#TOC li {
line-height: 16px;
}
table {
margin: 1em auto;
border-width: 1px;
border-color: #DDDDDD;
border-style: outset;
border-collapse: collapse;
}
table th {
border-width: 2px;
padding: 5px;
border-style: inset;
}
table td {
border-width: 1px;
border-style: inset;
line-height: 18px;
padding: 5px 5px;
}
table, table th, table td {
border-left-style: none;
border-right-style: none;
}
table thead, table tr.even {
background-color: #f7f7f7;
}
p {
margin: 0.5em 0;
}
blockquote {
background-color: #f6f6f6;
padding: 0.25em 0.75em;
}
hr {
border-style: solid;
border: none;
border-top: 1px solid #777;
margin: 28px 0;
}
dl {
margin-left: 0;
}
dl dd {
margin-bottom: 13px;
margin-left: 13px;
}
dl dt {
font-weight: bold;
}
ul {
margin-top: 0;
}
ul li {
list-style: circle outside;
}
ul ul {
margin-bottom: 0;
}
pre, code {
background-color: #f7f7f7;
border-radius: 3px;
color: #333;
white-space: pre-wrap; 
}
pre {
border-radius: 3px;
margin: 5px 0px 10px 0px;
padding: 10px;
}
pre:not([class]) {
background-color: #f7f7f7;
}
code {
font-family: Consolas, Monaco, 'Courier New', monospace;
font-size: 85%;
}
p > code, li > code {
padding: 2px 0px;
}
div.figure {
text-align: center;
}
img {
background-color: #FFFFFF;
padding: 2px;
border: 1px solid #DDDDDD;
border-radius: 3px;
border: 1px solid #CCCCCC;
margin: 0 5px;
}
h1 {
margin-top: 0;
font-size: 35px;
line-height: 40px;
}
h2 {
border-bottom: 4px solid #f7f7f7;
padding-top: 10px;
padding-bottom: 2px;
font-size: 145%;
}
h3 {
border-bottom: 2px solid #f7f7f7;
padding-top: 10px;
font-size: 120%;
}
h4 {
border-bottom: 1px solid #f7f7f7;
margin-left: 8px;
font-size: 105%;
}
h5, h6 {
border-bottom: 1px solid #ccc;
font-size: 105%;
}
a {
color: #0033dd;
text-decoration: none;
}
a:hover {
color: #6666ff; }
a:visited {
color: #800080; }
a:visited:hover {
color: #BB00BB; }
a[href^="http:"] {
text-decoration: underline; }
a[href^="https:"] {
text-decoration: underline; }

code > span.kw { color: #555; font-weight: bold; } 
code > span.dt { color: #902000; } 
code > span.dv { color: #40a070; } 
code > span.bn { color: #d14; } 
code > span.fl { color: #d14; } 
code > span.ch { color: #d14; } 
code > span.st { color: #d14; } 
code > span.co { color: #888888; font-style: italic; } 
code > span.ot { color: #007020; } 
code > span.al { color: #ff0000; font-weight: bold; } 
code > span.fu { color: #900; font-weight: bold; } 
code > span.er { color: #a61717; background-color: #e3d2d2; } 
</style>




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
<span id="cb1-2"><a href="#cb1-2" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Skipping install of &#39;ipcwQR&#39; from a github remote, the SHA1 (34689cfe) has not changed since last install.</span></span>
<span id="cb1-3"><a href="#cb1-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   Use `force = TRUE` to force installation</span></span></code></pre></div>
<p>or</p>
<div class="sourceCode" id="cb2"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a>base<span class="sc">::</span><span class="fu">require</span>(<span class="st">&quot;ipcwQR&quot;</span>)</span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: ipcwQR</span></span></code></pre></div>
<p>picrq() function has the following arguments:</p>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="co"># picrq(L=V,R=U,delta=delta,x=x,tau=tau)</span></span></code></pre></div>
<p>see the detailed description from help(picrq()).</p>
<p>We first simulate univariate partly interval-censored data with
normal random error, which is a similar simulation setting of Kim et
al. (2022+).</p>
<div class="sourceCode" id="cb4"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb4-1"><a href="#cb4-1" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(ipcwQR)</span>
<span id="cb4-2"><a href="#cb4-2" aria-hidden="true" tabindex="-1"></a><span class="fu">set.seed</span>(<span class="dv">316</span>)</span>
<span id="cb4-3"><a href="#cb4-3" aria-hidden="true" tabindex="-1"></a>n <span class="ot">=</span> <span class="dv">200</span></span>
<span id="cb4-4"><a href="#cb4-4" aria-hidden="true" tabindex="-1"></a>x1 <span class="ot">=</span> <span class="fu">runif</span>(n,<span class="sc">-</span><span class="dv">1</span>,<span class="dv">1</span>)</span>
<span id="cb4-5"><a href="#cb4-5" aria-hidden="true" tabindex="-1"></a>x2 <span class="ot">=</span> <span class="fu">rbinom</span>(n,<span class="dv">1</span>,<span class="fl">0.43</span>)</span>
<span id="cb4-6"><a href="#cb4-6" aria-hidden="true" tabindex="-1"></a>x <span class="ot">=</span> <span class="fu">cbind</span>(x1,x2)</span>
<span id="cb4-7"><a href="#cb4-7" aria-hidden="true" tabindex="-1"></a>T <span class="ot">=</span> <span class="dv">2</span> <span class="sc">+</span> x1 <span class="sc">+</span> x2 <span class="sc">+</span> <span class="fu">rnorm</span>(n)</span>
<span id="cb4-8"><a href="#cb4-8" aria-hidden="true" tabindex="-1"></a>U <span class="ot">=</span> (<span class="dv">1</span> <span class="sc">-</span> <span class="fl">0.25</span><span class="sc">*</span>x1)<span class="sc">*</span><span class="fu">runif</span>(n, <span class="sc">-</span><span class="dv">6</span>, <span class="dv">5</span>)</span>
<span id="cb4-9"><a href="#cb4-9" aria-hidden="true" tabindex="-1"></a>V <span class="ot">=</span> U <span class="sc">+</span> (<span class="dv">1</span> <span class="sc">-</span> <span class="fl">0.1</span><span class="sc">*</span>x2)<span class="sc">*</span><span class="fu">runif</span>(n, <span class="dv">6</span>, <span class="dv">20</span>)</span>
<span id="cb4-10"><a href="#cb4-10" aria-hidden="true" tabindex="-1"></a>U <span class="ot">=</span> <span class="fu">exp</span>(dplyr<span class="sc">::</span><span class="fu">case_when</span>(<span class="cn">TRUE</span> <span class="sc">~</span> T, T<span class="sc">&gt;</span>V <span class="sc">~</span> V, T<span class="sc">&lt;</span>U <span class="sc">~</span> <span class="sc">-</span><span class="cn">Inf</span>))</span>
<span id="cb4-11"><a href="#cb4-11" aria-hidden="true" tabindex="-1"></a>V <span class="ot">=</span> <span class="fu">exp</span>(dplyr<span class="sc">::</span><span class="fu">case_when</span>(<span class="cn">TRUE</span> <span class="sc">~</span> T, T<span class="sc">&gt;</span>V <span class="sc">~</span> <span class="cn">Inf</span>, T<span class="sc">&lt;</span>U <span class="sc">~</span> U))</span>
<span id="cb4-12"><a href="#cb4-12" aria-hidden="true" tabindex="-1"></a>delta <span class="ot">=</span> <span class="fu">ifelse</span>(U<span class="sc">==</span>V, <span class="dv">1</span>, <span class="dv">0</span>)</span>
<span id="cb4-13"><a href="#cb4-13" aria-hidden="true" tabindex="-1"></a>tau<span class="ot">=</span><span class="fl">0.3</span></span>
<span id="cb4-14"><a href="#cb4-14" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(<span class="at">L=</span>V,<span class="at">R=</span>U,<span class="at">delta=</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau)</span>
<span id="cb4-15"><a href="#cb4-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Warning in rq.fit.br(wx, wy, tau = tau, ...): Solution may be nonunique</span></span>
<span id="cb4-16"><a href="#cb4-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb4-17"><a href="#cb4-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.597118 0.018059      0 1.561723 1.632512</span></span>
<span id="cb4-18"><a href="#cb4-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     1.009903 0.030962      0 0.949217 1.070589</span></span>
<span id="cb4-19"><a href="#cb4-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     1.054783 0.061964      0 0.933334 1.176233</span></span>
<span id="cb4-20"><a href="#cb4-20" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(<span class="at">L=</span>V,<span class="at">R=</span>U,<span class="at">delta=</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">estimation =</span> <span class="st">&quot;dr&quot;</span>)</span>
<span id="cb4-21"><a href="#cb4-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Warning in rq.fit.br(wx, wy, tau = tau, ...): Solution may be nonunique</span></span>
<span id="cb4-22"><a href="#cb4-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se pvalue lower bd upper bd</span></span>
<span id="cb4-23"><a href="#cb4-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     1.597118 0.018059      0 1.561723 1.632513</span></span>
<span id="cb4-24"><a href="#cb4-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x1        0.3     1.009903 0.030962      0 0.949217 1.070589</span></span>
<span id="cb4-25"><a href="#cb4-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; x2        0.3     1.054785 0.061964      0 0.933335 1.176234</span></span></code></pre></div>
<p>We posit two estimating methods, Gehan and log-rank, which can be
conducted by specifying type = “gehan” and type = “logrank”,
respectively.</p>
<p>Next, we give illustrative example of the method for the multivariate
clustered data, by using phase 3 metastatic colorectal cancer clinical
trial. This dataset is available data(mCRC) in PICBayes package (Pan,
2021).</p>
<div class="sourceCode" id="cb5"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(ipcwQR)</span>
<span id="cb5-2"><a href="#cb5-2" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(tidyverse)</span>
<span id="cb5-3"><a href="#cb5-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──</span></span>
<span id="cb5-4"><a href="#cb5-4" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 </span></span>
<span id="cb5-5"><a href="#cb5-5" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ tibble  3.1.8      ✔ dplyr   1.0.10</span></span>
<span id="cb5-6"><a href="#cb5-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ tidyr   1.2.0      ✔ stringr 1.4.0 </span></span>
<span id="cb5-7"><a href="#cb5-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✔ readr   2.1.2      ✔ forcats 0.5.2 </span></span>
<span id="cb5-8"><a href="#cb5-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──</span></span>
<span id="cb5-9"><a href="#cb5-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✖ dplyr::filter() masks stats::filter()</span></span>
<span id="cb5-10"><a href="#cb5-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ✖ dplyr::lag()    masks stats::lag()</span></span>
<span id="cb5-11"><a href="#cb5-11" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(PICBayes)</span>
<span id="cb5-12"><a href="#cb5-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: coda</span></span>
<span id="cb5-13"><a href="#cb5-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: MCMCpack</span></span>
<span id="cb5-14"><a href="#cb5-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Loading required package: MASS</span></span>
<span id="cb5-15"><a href="#cb5-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb5-16"><a href="#cb5-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Attaching package: &#39;MASS&#39;</span></span>
<span id="cb5-17"><a href="#cb5-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb5-18"><a href="#cb5-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; The following object is masked from &#39;package:dplyr&#39;:</span></span>
<span id="cb5-19"><a href="#cb5-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb5-20"><a href="#cb5-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     select</span></span>
<span id="cb5-21"><a href="#cb5-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb5-22"><a href="#cb5-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb5-23"><a href="#cb5-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Markov Chain Monte Carlo Package (MCMCpack)</span></span>
<span id="cb5-24"><a href="#cb5-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Copyright (C) 2003-2022 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park</span></span>
<span id="cb5-25"><a href="#cb5-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb5-26"><a href="#cb5-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## Support provided by the U.S. National Science Foundation</span></span>
<span id="cb5-27"><a href="#cb5-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ## (Grants SES-0350646 and SES-0350613)</span></span>
<span id="cb5-28"><a href="#cb5-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; ##</span></span>
<span id="cb5-29"><a href="#cb5-29" aria-hidden="true" tabindex="-1"></a><span class="fu">data</span>(<span class="st">&quot;mCRC&quot;</span>)</span>
<span id="cb5-30"><a href="#cb5-30" aria-hidden="true" tabindex="-1"></a>d <span class="ot">=</span> <span class="fu">with</span>(<span class="fu">data.frame</span>(mCRC), <span class="fu">data.frame</span>(<span class="at">L =</span> <span class="fu">as.numeric</span>(L),</span>
<span id="cb5-31"><a href="#cb5-31" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">R =</span> <span class="fu">as.numeric</span>(R),</span>
<span id="cb5-32"><a href="#cb5-32" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">U =</span> <span class="fu">ifelse</span>(y<span class="sc">==</span><span class="dv">0</span>,R,L),</span>
<span id="cb5-33"><a href="#cb5-33" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">V =</span> <span class="fu">ifelse</span>(y<span class="sc">==</span><span class="dv">2</span>,L,R),</span>
<span id="cb5-34"><a href="#cb5-34" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Cluster weighted data</span></span>
<span id="cb5-35"><a href="#cb5-35" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">id=</span>(<span class="fu">rep</span>(<span class="fu">c</span>(<span class="fu">table</span>(SITE)),<span class="fu">c</span>(<span class="fu">table</span>(SITE)))),</span>
<span id="cb5-36"><a href="#cb5-36" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI.</span></span>
<span id="cb5-37"><a href="#cb5-37" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">x1=</span> <span class="fu">case_when</span>(TRT_C <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">0</span>, <span class="co">#Pan et al data</span></span>
<span id="cb5-38"><a href="#cb5-38" aria-hidden="true" tabindex="-1"></a>                                                    TRT_C <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">1</span>),</span>
<span id="cb5-39"><a href="#cb5-39" aria-hidden="true" tabindex="-1"></a>                                      <span class="co"># Tumor KRAS mutation status: 0 = wild-type, 1 = mutant.</span></span>
<span id="cb5-40"><a href="#cb5-40" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">x2=</span> <span class="fu">case_when</span>(KRAS_C <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">1</span>,</span>
<span id="cb5-41"><a href="#cb5-41" aria-hidden="true" tabindex="-1"></a>                                                    KRAS_C <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">0</span>),</span>
<span id="cb5-42"><a href="#cb5-42" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">site =</span> <span class="fu">as.numeric</span>(SITE),</span>
<span id="cb5-43"><a href="#cb5-43" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">y =</span> <span class="fu">as.numeric</span>(y),</span>
<span id="cb5-44"><a href="#cb5-44" aria-hidden="true" tabindex="-1"></a>                                      <span class="at">delta =</span> <span class="fu">case_when</span>(IC <span class="sc">==</span> <span class="dv">0</span> <span class="sc">~</span> <span class="dv">1</span>,</span>
<span id="cb5-45"><a href="#cb5-45" aria-hidden="true" tabindex="-1"></a>                                                        IC <span class="sc">==</span> <span class="dv">1</span> <span class="sc">~</span> <span class="dv">0</span>)</span>
<span id="cb5-46"><a href="#cb5-46" aria-hidden="true" tabindex="-1"></a>));</span>
<span id="cb5-47"><a href="#cb5-47" aria-hidden="true" tabindex="-1"></a>L<span class="ot">=</span>d<span class="sc">$</span>U;R<span class="ot">=</span>d<span class="sc">$</span>V; delta<span class="ot">=</span>d<span class="sc">$</span>delta</span>
<span id="cb5-48"><a href="#cb5-48" aria-hidden="true" tabindex="-1"></a>x <span class="ot">=</span> <span class="fu">cbind</span>(d<span class="sc">$</span>x1,d<span class="sc">$</span>x2); tau<span class="ot">=</span><span class="fl">0.3</span></span>
<span id="cb5-49"><a href="#cb5-49" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(d<span class="sc">$</span>U,d<span class="sc">$</span>V,d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau)</span>
<span id="cb5-50"><a href="#cb5-50" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb5-51"><a href="#cb5-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     3.309579 0.236915 0.000000  2.845225 3.773932</span></span>
<span id="cb5-52"><a href="#cb5-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.091733 0.131044 0.241957 -0.165112 0.348579</span></span>
<span id="cb5-53"><a href="#cb5-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.641756 0.187960 0.000320  0.273355 1.010157</span></span>
<span id="cb5-54"><a href="#cb5-54" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(d<span class="sc">$</span>U,d<span class="sc">$</span>V,d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau, <span class="at">estimation =</span> <span class="st">&quot;dr&quot;</span>)</span>
<span id="cb5-55"><a href="#cb5-55" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb5-56"><a href="#cb5-56" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     3.309579 0.236915 0.000000  2.845225 3.773932</span></span>
<span id="cb5-57"><a href="#cb5-57" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.091733 0.131044 0.241957 -0.165113 0.348579</span></span>
<span id="cb5-58"><a href="#cb5-58" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.641756 0.187960 0.000320  0.273355 1.010157</span></span>
<span id="cb5-59"><a href="#cb5-59" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(d<span class="sc">$</span>U,d<span class="sc">$</span>V,d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">wttype =</span> <span class="st">&quot;nonparam&quot;</span>,<span class="at">h=</span><span class="fl">0.9</span>)</span>
<span id="cb5-60"><a href="#cb5-60" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb5-61"><a href="#cb5-61" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     3.208453 0.220519 0.000000  2.776236 3.640670</span></span>
<span id="cb5-62"><a href="#cb5-62" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.192817 0.160906 0.115396 -0.122559 0.508192</span></span>
<span id="cb5-63"><a href="#cb5-63" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.641801 0.148979 0.000008  0.349803 0.933799</span></span>
<span id="cb5-64"><a href="#cb5-64" aria-hidden="true" tabindex="-1"></a><span class="fu">picrq</span>(d<span class="sc">$</span>U,d<span class="sc">$</span>V,d<span class="sc">$</span>delta,<span class="at">x=</span>x,<span class="at">tau=</span>tau,<span class="at">wttype =</span> <span class="st">&quot;nonparam&quot;</span>,<span class="at">id=</span>d<span class="sc">$</span>id,<span class="at">h=</span><span class="fl">0.9</span>)</span>
<span id="cb5-65"><a href="#cb5-65" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           tau coefficients       se   pvalue  lower bd upper bd</span></span>
<span id="cb5-66"><a href="#cb5-66" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Intercept 0.3     3.208297 0.220508 0.000000  2.776101 3.640494</span></span>
<span id="cb5-67"><a href="#cb5-67" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2         0.3     0.192902 0.160925 0.115322 -0.122512 0.508316</span></span>
<span id="cb5-68"><a href="#cb5-68" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3         0.3     0.641853 0.148961 0.000008  0.349888 0.933817</span></span></code></pre></div>
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



<!-- code folding -->


<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

</body>
</html>
