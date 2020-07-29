---
title: "The T-index: Measuring the reliability of accuracy estimates obtained from non-probability samples"
author: "Franz Waldner"
date: "15/07/2020"
output: html_document
---

# Introduction

This repository contains the code to compute the T index, an index that evalutates if non-probability samples provide reliable estimates of map accuracy.


# Context

In remote sensing, the term accuracy typically expresses the degree of correctness of a map.  Best practices in accuracy assessment have been widely researched  and include guidelines on how to select validation data using probability sampling designs. In practice, however,  probability samples data may be lacking and cross-validation based on non-probability samples is commonly used instead. This practice is risky because the resulting accuracy estimates can easily be mistaken for map accuracy. The following question arises: to what extent are accuracy estimates obtained from non-probability samples representative of map accuracy? This letter introduces the $T$ index to answer this question. Certain cross-validation designs (such as the common single-split or hold-out validation) provide representative accuracy estimates when hold-out sets are simple random samples of the map population. The $T$ index essentially measures the probability of a hold-out set of unknown sampling design to be a simple random sample. To that aim, we compare its spread in the feature space  against the spread of random unlabelled samples of the same size. The data spread is measured by a variant of Moran's $I$ autocorrelation index. Consistent interpretation of the $T$ index is proposed through the prism of significance testing, with $T$ values $<$ 0.05 indicating unreliable accuracy estimates. Its relevance and interpretation guidelines are also illustrated in a case study on crop-type mapping. Uptake of the $T$ index by the remote-sensing community will help inform about--sometimes caution against--the representativeness of accuracy estimates obtained by cross-validation, so that users can better decide whether a map is fit for their purpose or how its accuracy impacts their application. Subsequently, the $T$ index will build trust and improve the transparency of accuracy assessment in conditions which deviate from gold standards. 


<img src="./images/rationale.png" alt="T"
	title="Procedure to construct theTindex:  1) calculate the normalised Moran’sIindex of thelabelled hold-out set; 3) generate random unlabelled samples from the map population; 3) calculatethe normalised Moran’sIindex of all random unlabelled samples; 4) compute the probability of thelabelled set to belong to the empirical distribution of random unlabelled samples" width="350" />


# Structure of the repository

Functions to compute the the normalised Moran's $I$ index and the $T$ index are in the `R` folder.

Vignettes are under construction.


