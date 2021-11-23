---
title: "Eye size and visual range for bioluminescence"
author: "Katie Thomas"
date: "9/21/2021"
output:
  html_document:
    keep_md: true
    code_fold: hide
    theme: flatly
    toc: yes
    toc_float: yes
  pdf_document:
    toc: yes
---

<style type="text/css">

body{ /* Normal  */
      font-size: 17px;
  }
  
</style>




# Modeling of sighting distance for bioluminescent point source

Maximum sighting distances for point-source bioluminescence were modeled based on eye size across species. Results are imported and plotted here. 


```r
#import modeling results
modeling <- data.frame(read.csv("../Data/modeling.csv", header=TRUE, na.strings=c("", "NA", " ", stringsAsFactors = FALSE)))
```


```r
#reorder factor levels for figure legend (phylo order)
modeling$species <- factor(modeling$species, 
                                      levels = c("Deosergestes corniculum",
                                                 "Deosergestes henseni",
                                                 "Allosergestes pectinatus",
                                                 "Allosergestes sargassi",
                                                 "Sergestes atlanticus",
                                                 "Neosergestes edwardsii",
                                                 "Parasergestes vigilax",
                                                 "Parasergestes armatus",
                                                 "Eusergestes arcticus",
                                                 "Gardinerosergia splendens",
                                                 "Robustosergia regalis",
                                                 "Robustosergia robusta",
                                                 "Phorcosergia grandis",
                                                 "Sergia tenuiremis",
                                                 "Challengerosergia talismani",
                                                 "Challengerosergia hansjacobi"))

#make shape palette
shapes.sp <- c("Deosergestes corniculum" = 21,
              "Deosergestes henseni" = 22, 
              "Allosergestes pectinatus" = 23, 
              "Allosergestes sargassi" = 24,
              "Sergestes atlanticus" = 25, 
              "Neosergestes edwardsii"= 23,
              "Parasergestes vigilax" = 21, 
              "Parasergestes armatus" = 22,
              "Eusergestes arcticus" = 23, 
              "Gardinerosergia splendens" = 24, 
              "Robustosergia regalis" = 25, 
              "Robustosergia robusta" = 21, 
              "Phorcosergia grandis" = 22,
              "Sergia tenuiremis" = 23,
              "Challengerosergia talismani" = 24,
              "Challengerosergia hansjacobi" = 25)

#sergia/sergestes green purple pallette
cols.sp <- c(#Sergestes group
              "Deosergestes corniculum" = "#512E5F",
              "Deosergestes henseni" = "#633974",
              "Allosergestes pectinatus" = "#76448A",
              "Allosergestes sargassi" = "#884EA0",
              "Sergestes atlanticus" = "#9B59B6",
              "Neosergestes edwardsii" = "#AF7AC5",
              "Parasergestes vigilax" = "#C39BD3",
              "Parasergestes armatus" = "#D7BDE2",
              "Eusergestes arcticus" = "#EBDEF0",
             #Sergia group
              "Gardinerosergia splendens" = "#ABEBC6",
              "Robustosergia regalis" = "#A9DFBF",
              "Robustosergia robusta" = "#52BE80",
              "Phorcosergia grandis" = "#27AE60",
              "Sergia tenuiremis" = "#1E8449",
              "Challengerosergia talismani" = "#196F3D",
              "Challengerosergia hansjacobi" = "#145A32")
```

## Absolute sighting distance


```r
# OLS regression of absolute sighting distance v. body length
ols_abs <- lm(abs_sight_m ~ body_mm, data = modeling)

#model output
summary(ols_abs)
```

```{style="max-height: 300px;"}
## 
## Call:
## lm(formula = abs_sight_m ~ body_mm, data = modeling)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.5269 -0.2235  0.0165  0.1402  0.4943 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 0.340994   0.235249   1.449    0.169    
## body_mm     0.052807   0.006043   8.739 4.82e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3168 on 14 degrees of freedom
## Multiple R-squared:  0.8451,	Adjusted R-squared:  0.834 
## F-statistic: 76.37 on 1 and 14 DF,  p-value: 4.825e-07
```

```r
#save coefficient estimates
cc_olsabs <- as.data.frame(coefficients(ols_abs))

# Make plot for absolute sighting distances
plot_model_abs <- ggplot(modeling, aes(x = body_mm, y = abs_sight_m, color =  species, shape = species, fill = species)) + 
  geom_point(size = 2, alpha = 1) + 
  scale_shape_manual(values = shapes.sp, name = "Species") + 
  scale_color_manual(values = cols.sp, name = "Species") +
  scale_fill_manual(values = cols.sp, name = "Species") +
  xlab("Body length (mm)") + 
  ylab("Sighting distance (m)") + 
  xlim(c(0,60)) +
  ylim(c(0,4)) +
  theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(face = "italic")) +
   geom_abline(data = cc_olsabs, aes(intercept = cc_olsabs[1,1], slope = cc_olsabs[2,1]), linetype = "dashed")

# Interactive plot
ggplotly(plot_model_abs)
```

```{=html}
<div id="htmlwidget-5f6fea9dc6039b3d41ef" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-5f6fea9dc6039b3d41ef">{"x":{"data":[{"x":[48.59],"y":[2.38],"text":"body_mm: 48.59<br />abs_sight_m: 2.38<br />species: Deosergestes corniculum<br />species: Deosergestes corniculum<br />species: Deosergestes corniculum","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(81,46,95,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(81,46,95,1)"}},"hoveron":"points","name":"Deosergestes corniculum","legendgroup":"Deosergestes corniculum","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[40.43],"y":[2.17],"text":"body_mm: 40.43<br />abs_sight_m: 2.17<br />species: Deosergestes henseni<br />species: Deosergestes henseni<br />species: Deosergestes henseni","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(99,57,116,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(99,57,116,1)"}},"hoveron":"points","name":"Deosergestes henseni","legendgroup":"Deosergestes henseni","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[19.98],"y":[1.11],"text":"body_mm: 19.98<br />abs_sight_m: 1.11<br />species: Allosergestes pectinatus<br />species: Allosergestes pectinatus<br />species: Allosergestes pectinatus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(118,68,138,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(118,68,138,1)"}},"hoveron":"points","name":"Allosergestes pectinatus","legendgroup":"Allosergestes pectinatus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[23.82],"y":[1.54],"text":"body_mm: 23.82<br />abs_sight_m: 1.54<br />species: Allosergestes sargassi<br />species: Allosergestes sargassi<br />species: Allosergestes sargassi","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(136,78,160,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(136,78,160,1)"}},"hoveron":"points","name":"Allosergestes sargassi","legendgroup":"Allosergestes sargassi","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[25.44],"y":[1.76],"text":"body_mm: 25.44<br />abs_sight_m: 1.76<br />species: Sergestes atlanticus<br />species: Sergestes atlanticus<br />species: Sergestes atlanticus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(155,89,182,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(155,89,182,1)"}},"hoveron":"points","name":"Sergestes atlanticus","legendgroup":"Sergestes atlanticus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[16.89],"y":[1.33],"text":"body_mm: 16.89<br />abs_sight_m: 1.33<br />species: Neosergestes edwardsii<br />species: Neosergestes edwardsii<br />species: Neosergestes edwardsii","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(175,122,197,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(175,122,197,1)"}},"hoveron":"points","name":"Neosergestes edwardsii","legendgroup":"Neosergestes edwardsii","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[19.12],"y":[1.33],"text":"body_mm: 19.12<br />abs_sight_m: 1.33<br />species: Parasergestes vigilax<br />species: Parasergestes vigilax<br />species: Parasergestes vigilax","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(195,155,211,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(195,155,211,1)"}},"hoveron":"points","name":"Parasergestes vigilax","legendgroup":"Parasergestes vigilax","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[30.64],"y":[1.54],"text":"body_mm: 30.64<br />abs_sight_m: 1.54<br />species: Parasergestes armatus<br />species: Parasergestes armatus<br />species: Parasergestes armatus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(215,189,226,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(215,189,226,1)"}},"hoveron":"points","name":"Parasergestes armatus","legendgroup":"Parasergestes armatus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[54],"y":[2.99],"text":"body_mm: 54.00<br />abs_sight_m: 2.99<br />species: Eusergestes arcticus<br />species: Eusergestes arcticus<br />species: Eusergestes arcticus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(235,222,240,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(235,222,240,1)"}},"hoveron":"points","name":"Eusergestes arcticus","legendgroup":"Eusergestes arcticus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[33.04],"y":[2.58],"text":"body_mm: 33.04<br />abs_sight_m: 2.58<br />species: Gardinerosergia splendens<br />species: Gardinerosergia splendens<br />species: Gardinerosergia splendens","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(171,235,198,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(171,235,198,1)"}},"hoveron":"points","name":"Gardinerosergia splendens","legendgroup":"Gardinerosergia splendens","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[45.06],"y":[2.99],"text":"body_mm: 45.06<br />abs_sight_m: 2.99<br />species: Robustosergia regalis<br />species: Robustosergia regalis<br />species: Robustosergia regalis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,223,191,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(169,223,191,1)"}},"hoveron":"points","name":"Robustosergia regalis","legendgroup":"Robustosergia regalis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[57.63],"y":[3.77],"text":"body_mm: 57.63<br />abs_sight_m: 3.77<br />species: Robustosergia robusta<br />species: Robustosergia robusta<br />species: Robustosergia robusta","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(82,190,128,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(82,190,128,1)"}},"hoveron":"points","name":"Robustosergia robusta","legendgroup":"Robustosergia robusta","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[56.07],"y":[3.19],"text":"body_mm: 56.07<br />abs_sight_m: 3.19<br />species: Phorcosergia grandis<br />species: Phorcosergia grandis<br />species: Phorcosergia grandis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(39,174,96,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(39,174,96,1)"}},"hoveron":"points","name":"Phorcosergia grandis","legendgroup":"Phorcosergia grandis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[45.36],"y":[2.79],"text":"body_mm: 45.36<br />abs_sight_m: 2.79<br />species: Sergia tenuiremis<br />species: Sergia tenuiremis<br />species: Sergia tenuiremis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(30,132,73,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(30,132,73,1)"}},"hoveron":"points","name":"Sergia tenuiremis","legendgroup":"Sergia tenuiremis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[37.14],"y":[2.38],"text":"body_mm: 37.14<br />abs_sight_m: 2.38<br />species: Challengerosergia talismani<br />species: Challengerosergia talismani<br />species: Challengerosergia talismani","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(25,111,61,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(25,111,61,1)"}},"hoveron":"points","name":"Challengerosergia talismani","legendgroup":"Challengerosergia talismani","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[33.34],"y":[2.58],"text":"body_mm: 33.34<br />abs_sight_m: 2.58<br />species: Challengerosergia hansjacobi<br />species: Challengerosergia hansjacobi<br />species: Challengerosergia hansjacobi","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(20,90,50,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(20,90,50,1)"}},"hoveron":"points","name":"Challengerosergia hansjacobi","legendgroup":"Challengerosergia hansjacobi","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-3,-3,63,63],"y":[0.182572262634884,0.182572262634884,3.66785110694591,3.66785110694591],"text":"cc_olsabs[1, 1]: 0.340994<br />cc_olsabs[2, 1]: 0.05280726","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":27.689497716895,"r":7.30593607305936,"b":41.6438356164384,"l":31.4155251141553},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-3,63],"tickmode":"array","ticktext":["0","20","40","60"],"tickvals":[0,20,40,60],"categoryorder":"array","categoryarray":["0","20","40","60"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Body length (mm)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.2,4.2],"tickmode":"array","ticktext":["0","1","2","3","4"],"tickvals":[0,1,2,3,4],"categoryorder":"array","categoryarray":["0","1","2","3","4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"Sighting distance (m)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"title":{"text":"Species","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"21b91fc336d2":{"x":{},"y":{},"colour":{},"shape":{},"fill":{},"type":"scatter"},"21b949442548":{"intercept":{},"slope":{}}},"cur_data":"21b91fc336d2","visdat":{"21b91fc336d2":["function (y) ","x"],"21b949442548":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

## Sighting distance relative to body length


```r
# OLS regression of relative sighting distance v. body length
ols_rel <- lm(rel_sight ~ body_mm, data = modeling)

#model output
summary(ols_rel)
```

```{style="max-height: 300px;"}
## 
## Call:
## lm(formula = rel_sight ~ body_mm, data = modeling)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -14.8898  -4.4492   0.5978   5.9812  13.6333 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  74.4142     6.6957  11.114 2.49e-08 ***
## body_mm      -0.2971     0.1720  -1.728    0.106    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9.016 on 14 degrees of freedom
## Multiple R-squared:  0.1757,	Adjusted R-squared:  0.1169 
## F-statistic: 2.985 on 1 and 14 DF,  p-value: 0.106
```

```r
#save coefficient estimates
cc_olsrel <- as.data.frame(coefficients(ols_rel))

# Make plot for relative sighting distances
plot_model_rel <- ggplot(modeling, aes(x = body_mm, y = rel_sight, color =  species, shape = species, fill = species)) + 
  geom_point(size = 2, alpha = 1) + 
  scale_shape_manual(values = shapes.sp, name = "Species") + 
  scale_color_manual(values = cols.sp, name = "Species") +
  scale_fill_manual(values = cols.sp, name = "Species") +
  xlab("Body length (mm)") + 
  ylab("Relative sighting distance") + 
  xlim(c(0,60)) +
  ylim(c(0,90)) +
  theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(face = "italic")) +
   geom_abline(data = cc_olsrel, aes(intercept = cc_olsrel[1,1], slope = cc_olsrel[2,1]), linetype = "dashed")

# Interactive plot
ggplotly(plot_model_rel)
```

```{=html}
<div id="htmlwidget-7893d005e3429edcf16a" style="width:768px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-7893d005e3429edcf16a">{"x":{"data":[{"x":[48.59],"y":[48.99],"text":"body_mm: 48.59<br />rel_sight: 48.99<br />species: Deosergestes corniculum<br />species: Deosergestes corniculum<br />species: Deosergestes corniculum","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(81,46,95,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(81,46,95,1)"}},"hoveron":"points","name":"Deosergestes corniculum","legendgroup":"Deosergestes corniculum","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[40.43],"y":[53.78],"text":"body_mm: 40.43<br />rel_sight: 53.78<br />species: Deosergestes henseni<br />species: Deosergestes henseni<br />species: Deosergestes henseni","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(99,57,116,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(99,57,116,1)"}},"hoveron":"points","name":"Deosergestes henseni","legendgroup":"Deosergestes henseni","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[19.98],"y":[55.78],"text":"body_mm: 19.98<br />rel_sight: 55.78<br />species: Allosergestes pectinatus<br />species: Allosergestes pectinatus<br />species: Allosergestes pectinatus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(118,68,138,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(118,68,138,1)"}},"hoveron":"points","name":"Allosergestes pectinatus","legendgroup":"Allosergestes pectinatus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[23.82],"y":[64.85],"text":"body_mm: 23.82<br />rel_sight: 64.85<br />species: Allosergestes sargassi<br />species: Allosergestes sargassi<br />species: Allosergestes sargassi","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(136,78,160,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(136,78,160,1)"}},"hoveron":"points","name":"Allosergestes sargassi","legendgroup":"Allosergestes sargassi","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[25.44],"y":[69.05],"text":"body_mm: 25.44<br />rel_sight: 69.05<br />species: Sergestes atlanticus<br />species: Sergestes atlanticus<br />species: Sergestes atlanticus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(155,89,182,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(155,89,182,1)"}},"hoveron":"points","name":"Sergestes atlanticus","legendgroup":"Sergestes atlanticus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[16.89],"y":[78.78],"text":"body_mm: 16.89<br />rel_sight: 78.78<br />species: Neosergestes edwardsii<br />species: Neosergestes edwardsii<br />species: Neosergestes edwardsii","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(175,122,197,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(175,122,197,1)"}},"hoveron":"points","name":"Neosergestes edwardsii","legendgroup":"Neosergestes edwardsii","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[19.12],"y":[69.59],"text":"body_mm: 19.12<br />rel_sight: 69.59<br />species: Parasergestes vigilax<br />species: Parasergestes vigilax<br />species: Parasergestes vigilax","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(195,155,211,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(195,155,211,1)"}},"hoveron":"points","name":"Parasergestes vigilax","legendgroup":"Parasergestes vigilax","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[30.64],"y":[50.42],"text":"body_mm: 30.64<br />rel_sight: 50.42<br />species: Parasergestes armatus<br />species: Parasergestes armatus<br />species: Parasergestes armatus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(215,189,226,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(215,189,226,1)"}},"hoveron":"points","name":"Parasergestes armatus","legendgroup":"Parasergestes armatus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[54],"y":[55.31],"text":"body_mm: 54.00<br />rel_sight: 55.31<br />species: Eusergestes arcticus<br />species: Eusergestes arcticus<br />species: Eusergestes arcticus","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(235,222,240,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(235,222,240,1)"}},"hoveron":"points","name":"Eusergestes arcticus","legendgroup":"Eusergestes arcticus","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[33.04],"y":[78.23],"text":"body_mm: 33.04<br />rel_sight: 78.23<br />species: Gardinerosergia splendens<br />species: Gardinerosergia splendens<br />species: Gardinerosergia splendens","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(171,235,198,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(171,235,198,1)"}},"hoveron":"points","name":"Gardinerosergia splendens","legendgroup":"Gardinerosergia splendens","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[45.06],"y":[66.29],"text":"body_mm: 45.06<br />rel_sight: 66.29<br />species: Robustosergia regalis<br />species: Robustosergia regalis<br />species: Robustosergia regalis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(169,223,191,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(169,223,191,1)"}},"hoveron":"points","name":"Robustosergia regalis","legendgroup":"Robustosergia regalis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[57.63],"y":[65.42],"text":"body_mm: 57.63<br />rel_sight: 65.42<br />species: Robustosergia robusta<br />species: Robustosergia robusta<br />species: Robustosergia robusta","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(82,190,128,1)","opacity":1,"size":7.55905511811024,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(82,190,128,1)"}},"hoveron":"points","name":"Robustosergia robusta","legendgroup":"Robustosergia robusta","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[56.07],"y":[56.81],"text":"body_mm: 56.07<br />rel_sight: 56.81<br />species: Phorcosergia grandis<br />species: Phorcosergia grandis<br />species: Phorcosergia grandis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(39,174,96,1)","opacity":1,"size":7.55905511811024,"symbol":"square","line":{"width":1.88976377952756,"color":"rgba(39,174,96,1)"}},"hoveron":"points","name":"Phorcosergia grandis","legendgroup":"Phorcosergia grandis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[45.36],"y":[61.43],"text":"body_mm: 45.36<br />rel_sight: 61.43<br />species: Sergia tenuiremis<br />species: Sergia tenuiremis<br />species: Sergia tenuiremis","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(30,132,73,1)","opacity":1,"size":7.55905511811024,"symbol":"diamond","line":{"width":1.88976377952756,"color":"rgba(30,132,73,1)"}},"hoveron":"points","name":"Sergia tenuiremis","legendgroup":"Sergia tenuiremis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[37.14],"y":[64.08],"text":"body_mm: 37.14<br />rel_sight: 64.08<br />species: Challengerosergia talismani<br />species: Challengerosergia talismani<br />species: Challengerosergia talismani","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(25,111,61,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-up","line":{"width":1.88976377952756,"color":"rgba(25,111,61,1)"}},"hoveron":"points","name":"Challengerosergia talismani","legendgroup":"Challengerosergia talismani","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[33.34],"y":[77.53],"text":"body_mm: 33.34<br />rel_sight: 77.53<br />species: Challengerosergia hansjacobi<br />species: Challengerosergia hansjacobi<br />species: Challengerosergia hansjacobi","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(20,90,50,1)","opacity":1,"size":7.55905511811024,"symbol":"triangle-down","line":{"width":1.88976377952756,"color":"rgba(20,90,50,1)"}},"hoveron":"points","name":"Challengerosergia hansjacobi","legendgroup":"Challengerosergia hansjacobi","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-3,-3,63,63],"y":[75.3056148625494,75.3056148625494,55.6944111830093,55.6944111830093],"text":"cc_olsrel[1, 1]: 74.4142<br />cc_olsrel[2, 1]: -0.2971394","type":"scatter","mode":"lines","line":{"width":1.88976377952756,"color":"rgba(0,0,0,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":27.689497716895,"r":7.30593607305936,"b":41.6438356164384,"l":37.2602739726027},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-3,63],"tickmode":"array","ticktext":["0","20","40","60"],"tickvals":[0,20,40,60],"categoryorder":"array","categoryarray":["0","20","40","60"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Body length (mm)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-4.5,94.5],"tickmode":"array","ticktext":["0","25","50","75"],"tickvals":[0,25,50,75],"categoryorder":"array","categoryarray":["0","25","50","75"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"Relative sighting distance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"title":{"text":"Species","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"21b955120f7e":{"x":{},"y":{},"colour":{},"shape":{},"fill":{},"type":"scatter"},"21b914b744d4":{"intercept":{},"slope":{}}},"cur_data":"21b955120f7e","visdat":{"21b955120f7e":["function (y) ","x"],"21b914b744d4":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
#name panels
fig.a <- plot_model_abs + theme(legend.position = "none") 
fig.b <- plot_model_rel + theme(legend.position = "none") 

#arrange plots in panels
plots <- plot_grid(fig.a, fig.b,
           align = 'vh', 
           axis = 'lb',
           labels = c("A", "B"), #panel labels for figure
           hjust = -1, #adjustment for panel labels
           nrow = 1) #number of rows in grids

# extract legend from Rana temporaria figure
leg <- get_legend(fig.a + theme(legend.position="right"))

#export figure
pdf("../Figures/Fig-6.pdf", width = 12, height = 4.5)
plot_grid(plots, leg, ncol = 2, rel_widths = c(1, .3))
dev.off()
```

```{style="max-height: 300px;"}
## quartz_off_screen 
##                 2
```

