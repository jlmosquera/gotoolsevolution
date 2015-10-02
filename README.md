# Study of Evolution and Clustering of GO Tools Classified in SerbGO Database

__Contents__

1. Motivation
2. Data
3. Statistical Methods
4. Licensing
5. NOTE
6. Authors

## 1. Motivation

In recent years, many tools have been developed to assist with the analysis of experimental results based on the _Gene Ontology_ (GO) ([1]). Some efforts to classify these tools based on the type of method applied for the enrichment of GO terms have been performed ([2]). But in fact, when a researcher is looking for software to perform an enrichment analysis, it is highly possible that it is lost, or at least that he/she does not find the most suitable tool for his/her needs, due to the large number of existing GO tools, even having classified the GO tools based on the type of enrichment that they carry out. In order to shed light on this issue, a web tool called [SerbGO](http://estbioinfo.stat.ub.es/apli/serbgov131/index.php) was developed to facilitate both searching for and comparing GO tools ([3]). Based on a _Standard Functionalities Set_, SerbGO database stores the classification of a long list of GO tools for enrichment analysis. 

The speed of development of new methods and tools for the enrichment analysis, as well as the improvement of existing applications, by the scientific community, is really considerable. In this regard, has any comprehensive monitoring of tools been not conducted to see how their capabilities evolve. In order to shed light on this issue, periodically, but not regularly, the SerbGO database is reviewed. This fact led us to observe a certain degree of evolution in tools classified in the SerbGO database. For this reason we decided to conduct a statistical study, based on the monitoring of all the tools included in the first version of SerbGO, in order to discuss the evolution of the functionalities, and to observe if some form of clustering of GO tools exists. 

## 2. Data

Data analysis is based on the number of standard functionalities that 26 tools had available in 3 different years (2005, 2007 and 2009). These tools are the list of the original GO tools stored in the first SerbGO database. For list of GO tools, six of the tables (type, species, data, annotation, statistics and outputs) stored in SerbGO database, corresponding to each
year, were downloaded.

After downloading information from SerbGO, an homogenization process of raw data has been performed in order to reduce redundancies.

The homogenization of these tables resulted in a reduction of the number of functionalities. From the original 205 standard functionalities, 178 functionalities have been selected and pre-processed to be analyzed.

For each year, the seven homogenized tables were merged into one unique matrix of large data. That is, data analysis has been based on three binary matrices such that each matrix describes the capabilities of the GO tools under study through the functionalities selected for a specific year (i.e. 2005, 2007 or 2009).

## 3. Statistical Methods

The statistical analysis to study the evolution and clustering of functionalities of GO tools has consisted of: 

* _Descriptive Statistics_: to provide basic summaries about samples and observations. It may be used to describe relationships between pairs of variables. 
* _Inferential Analysis_: to test hypotheses such as whether or not differences exist between frequencies of functionalities from different years.
* _Multivariate analyses_: to explore the behavior of GO tools according to their capabilities over time. That is, multivariate
methods have been applied to study how the different programs are grouped according to their capabilities throughout the years.

The data analysis has been performed with the statistical software R ([4]) supported by some extra packages, and explicitly programmed functions.

Raw data and functions used to perform the statistical analysis are available in this repository. 

## 4. Licensing

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">
 <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" />
</a>
<br />
<span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Study of Evolution and Clustering of GO Tools Classified in SerbGO Database</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Jose Luis Mosquera</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/jlmosquera/gotoolsevolution" rel="dct:source">https://github.com/jlmosquera/gotoolsevolution</a>.

## 5. NOTE

Data and R code provided in this repository are a part of the PhD thesis called *_Methods and Models for the Analysis of Biological Significance Based on High-Throughput Data_* at the [University of Barcelona](http://www.ub.edu/web/ub/en/index.html?).

## 6. Authors

* Author: Jose Luis Mosquera, PhD (jlmosquera@gmail.com) - Department of Statistics. Universtiy of Barcelona.
          Currently at the Bioinformatics/Biostatistics Uniat - IRB Barcelona.
* Advisor: Alex Sánchez, PhD - Department of Statistics. University of Barcelona.

## 7. References

1. Sorin Draghici, Purvesh Khatri, Rui P. Martins, G. Charles Ostermeier, and Stephen A. Krawetz. _Global Functional Profiling of Gene Gxpression._ Genomics, 81(2):98–104, 2003.
2. Da Wei Huang, Brad T. Sherman, and Richard A. Lempicki. _Bioinformatics Enrichment Tools: Paths Toward the Comprehensive Functional Analysis of Large Gene Lists. Nucleic Acids Research, 37(1):1–13, 2009.
3. Jose Luis Mosquera and Alex Sànchez-Pla. _SerbGO: Searching for the Best GO Tool._ Nucleic Acids Res., 36(Web Server Issue):W368–371, 2008.
4. R Core Team. _R: A Language and Environment for Statistical Computing._ R Foundation for Statistical Computing, Vienna, Austria, 2013.


