\name{samr-internal}
\title{Internal samr functions}
\alias{ttest.func}
\alias{wilcoxon.func}
\alias{onesample.ttest.func}
\alias{paired.ttest.func}
\alias{cox.func}
\alias{multiclass.func}
\alias{quantitative.func}
\alias{patterndiscovery.func}
\alias{timearea.func}
\alias{sumlengths}
\alias{detec.slab}
\alias{detec.horiz}
\alias{localfdr}
\alias{predictlocalfdr}
\alias{qvalue.func}
\alias{foldchange.twoclass}
\alias{foldchange.paired}
\alias{est.s0}
\alias{coxscor}
\alias{coxvar}
\alias{coxstuff}
\alias{permute}
\alias{permute.rows}
\alias{sample.perms}
\alias{insert.value}
\alias{integer.base.b}
\alias{varr}
\alias{mysvd}
\alias{getperms}
\alias{compute.block.perms}
\alias{parse.block.labels.for.2classes}
\alias{parse.time.labels.and.summarize.data}
\alias{check.format}
\alias{samr.const.quantitative.response}
\alias{samr.const.twoclass.unpaired.response}
\alias{samr.const.survival.response}
\alias{samr.const.multiclass.response}
\alias{samr.const.oneclass.response}
\alias{samr.const.twoclass.paired.response}
\alias{samr.const.twoclass.unpaired.timecourse.response}
\alias{samr.const.twoclass.paired.timecourse.response}
\alias{samr.const.oneclass.timecourse.response}
\alias{samr.const.patterndiscovery.response}
\alias{samr.xl.const.response.types.table}
\alias{samr.xl.build.data}
\alias{samr.xl.compute.plot.xy}
\alias{samr.xl.parse.survival}
\alias{samr.xl.parse.patterndiscovery}
\alias{samr.xl.error.trace}
\alias{samr.xl.impute.data}
\description{Internal samr functions}
\usage{
ttest.func(x, y, s0 = 0, sd=NULL)
wilcoxon.func(x, y, s0 = 0)
onesample.ttest.func(x, y, s0 = 0, sd=NULL)
paired.ttest.func(x, y, s0 = 0, sd=NULL)
cox.func(x, y, censoring.status, s0 = 0)
multiclass.func(x, y, s0 = 0)
quantitative.func(x, y, s0 = 0)
patterndiscovery.func(x,s0=0, eigengene.number=1)
timearea.func(x, y, s0 = 0)
sumlengths(aa)
detec.slab(samr.obj, del, min.foldchange)
detec.horiz(samr.obj, cutlow, cutup, min.foldchange)
localfdr(samr.obj, min.foldchange, perc = 0.01, df = 10)
predictlocalfdr(smooth.object, d)
qvalue.func(samr.obj, sig, delta.table)
foldchange.twoclass(x, y, logged2)
foldchange.paired(x, y, logged2)
est.s0(tt, sd, s0.perc = seq(0, 1, by = 0.05))
coxscor(x, y, ic, offset = rep(0, length(y)))
coxvar(x, y, ic, offset = rep(0, length(y)), coxstuff.obj = NULL)
coxstuff(x, y, ic, offset = rep(0, length(y)))
permute(elem)
sample.perms(elem, nperms)
permute.rows(x)
insert.value(vec, newval, pos)
integer.base.b(x, b = 2)
varr(x, meanx = NULL)
mysvd(x, n.components=NULL)
getperms(y, nperms)
compute.block.perms(y, blocky, nperms)
parse.block.labels.for.2classes(y)
parse.time.labels.and.summarize.data(x, y, resp.type,time.summary.type)

check.format(y, resp.type, censoring.status = NULL)
samr.const.quantitative.response
samr.const.twoclass.unpaired.response
samr.const.survival.response
samr.const.multiclass.response
samr.const.oneclass.response
samr.const.twoclass.paired.response
samr.const.twoclass.unpaired.timecourse.response
samr.const.twoclass.paired.timecourse.response
samr.const.oneclass.timecourse.response
samr.const.patterndiscovery.response
samr.xl.const.response.types.table
samr.xl.build.data(x, y, geneid, genenames, logged2, resp.type)
samr.xl.compute.plot.xy(samr.obj, del, min.foldchange = 0, plot = FALSE)
samr.xl.parse.survival(y)
samr.xl.parse.patterndiscovery(y)
samr.xl.error.trace()
samr.xl.impute.data(data.obj)
}
\author{Balasubrimanian Narasimhan and Robert Tibshirani}

\details{
  These are not to be called by the user.
}
\keyword{internal}

