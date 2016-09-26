# `platexpress` - Microbial Growth & Gene Expression


The platexpress package provides a quick & easy interface to microbial
growth & gene expression data as measured in parallel growth platforms
such as microplate readers. It allows for quick summarization over
replicates and quick comparison between experiments.

(TODO) A few data conversion routines allow to interface other R
packages for analysis of microbial growth such as `grofit` and
`growthcurver` and quickly display their results within `platexpress`.


## A Typical Workflow in platexpress
### 1) parse the plate layout map

```R
plate.file <- system.file("extdata", "AP12_layout.csv", package = "platexpress")
plate <- readPlateMap(file=plate.file, blank.id="blank",fsep="\n", fields=c("strain","samples"))
```

### 2) parse the data, exported from platereader software

```R
data.file <- system.file("extdata", "AP12.csv", package = "platexpress")
raw <- readPlateData(file=data.file, type="Synergy", data.ids=c("600","YFP_50:500,535"), dec=",")
```

### 3) inspect the raw data

```R
vp <- viewPlate(raw)
```

### 4) Note that there is no growth in A9, so let's skip it

```R
raw <- skipWells(raw, skip="A9")
```

### 5) Now correct for blank well measurements, and view only present
### rows/cols

```R
data <- correctBlanks(data=raw, plate=plate)
vp <- viewPlate(data, rows=c("A","B","C"),cols=1:9)
```

### 6) group replicates and view summarized growth/exprssion curves

```R
groups <- getGroups(plate, by=c("strain","samples"))
vg <- viewGroups(data,groups=groups,lwd.orig=0.5,nrow=3)
```
