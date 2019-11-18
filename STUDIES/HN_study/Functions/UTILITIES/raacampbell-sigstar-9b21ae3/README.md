# sigstar

Lines and asterisks indicating significant differences between two groups on a plot are commonly used in the life and social sciences. 
To my knowledge, no MATLAB function for adding these is openly available.
`sigstar` makes it easy to add lines and significance asterisks joining one or more pairs of groups on bar charts, box plots, and even line plots. The user simply supplies the identities of the two groups and the p-value 
(which the user has calculated using an appropriate test). 

<img src="https://github.com/raacampbell/sigstar/blob/gh-pages/images/ss.png" />

### Usage
Group identity is defined as x-axis locations or, optionally, group names (if the x-axis labels are strings). 
`sigstar` converts the supplied p-values to the appropriate number of asterisks and plots these over the lines that link the pair of groups.
`sigstar` attempts to intelligently place lines and asterisks so that they do not overlie existing plot elements. 
By default, bars with longer horizontal extents are plotted above shorter bars. 
This is the convention and it looks neater. 
The user has control over the order in which significance bars are added. 
Modifying the order provides control over the vertical position of the bars. 
This is important for obtaining a neat result in a cluttered plot.

### More
See examples in "`help sigstar`" for details.
The handles of the added plot elements are returned by default, providing the user with fine control of the plot's appearance. 
This is important since it's difficult to provision for all possible usage scenarios (see `demo_sigstar`).
The function should produce publication quality results, but you may need to play with the figure size and asterisk font size.

### Correct usage of significance stars
Significance indicators, like significance stars, are often missued. 
Keep the following guidelines in mind when using this function:
* Only use the stars to highlight important results. Do not clutter your figures with irrelevent statistical tests. 
* If your data can be fitted with a line or curve then do this and report and the coefficients and their associated confidence intervals. Do not perform multiple significance tests between different points along the curve. This is almost always bad practice. 
* Barplots hide data: avoid them (even with error bars). Instead use historgrams, violin plots, box plots, or [notBoxPlots](https://github.com/raacampbell/notBoxPlot).
