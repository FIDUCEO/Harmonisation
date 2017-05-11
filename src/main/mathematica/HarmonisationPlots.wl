(* ::Package:: *)

(* ::Title:: *)
(*A Wolfram Language Package to Create Harmonisation Plots with Mathematica*)


(* ::Text:: *)
(*This work was undertaken within the project "Fidelity and Uncertainty in Climate data records from Earth Observation" (FIDUCEO). FIDUCEO has received funding from the European Union's H2020 Research and Innovation programme, under Grant Agreement 638822.*)


(* ::Section:: *)
(*Plot functions*)


(* ::Input::Initialization:: *)
plotCorrelationMatrix[data_,OptionsPattern[{ColorFunction->"TemperatureMap",Frame->True,FrameLabel->{"Sensor","Sensor"},ParameterGrouping->1,PlotTheme->"Detailed"}]]:=Module[{labels,matrix,parameterGrouping,frameTicks},matrix=data[[1]];labels=data[[2]];parameterGrouping=OptionValue[ParameterGrouping];frameTicks=Table[{(parameterGrouping+1)/2+i-1,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {1,}"]->""]),45 Degree]},{i,1,Length[labels],parameterGrouping}];MatrixPlot[matrix,Mesh->{Length[matrix]/parameterGrouping-1,Length[matrix]/parameterGrouping-1},ColorFunction->OptionValue[ColorFunction],Frame->OptionValue[Frame],FrameLabel->OptionValue[FrameLabel],FrameTicks->{frameTicks,frameTicks},PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotCovarianceMatrix[data_,OptionsPattern[{ColorFunction->"TemperatureMap",Frame->True,FrameLabel->{"Sensor","Sensor"},ParameterGrouping->1,PlotTheme->"Detailed"}]]:=Module[{labels,matrix,parameterGrouping,frameTicks},matrix=data[[1]];labels=data[[2]];parameterGrouping=OptionValue[ParameterGrouping];frameTicks=Table[{(parameterGrouping+1)/2+i-1,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {1,}"]->""]),45 Degree]},{i,1,Length[labels],parameterGrouping}];MatrixPlot[matrix,Mesh->{Length[matrix]/parameterGrouping-1,Length[matrix]/parameterGrouping-1},ColorFunction->OptionValue[ColorFunction],Frame->OptionValue[Frame],FrameLabel->OptionValue[FrameLabel],FrameTicks->{frameTicks,frameTicks},PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotMeanResiduals[data_,OptionsPattern[{Exaggerate->100,PlotRange->Automatic,PlotTheme->{"Detailed","OpenMarkersThick","VibrantColor"}}]]:=Module[{exaggerate,plotRange,mean,stdev,range,labels,frameTicks},mean=data[[1]];exaggerate=OptionValue[Exaggerate];
plotRange=OptionValue[PlotRange];stdev=data[[2]];range=Max[Abs[{mean-exaggerate stdev,mean+exaggerate stdev}]];labels=data[[3]];frameTicks=Table[{i,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {2,}"]->""]),90 Degree]},{i,Length[labels]}];ListPlot[{mean,mean-exaggerate stdev,mean+exaggerate stdev},Filling->{1->{2},1->{3}},FillingStyle->Thick,FrameLabel->{"Sensors","Mean K residual ([radiance])"},FrameTicks->{{Automatic,Automatic},{frameTicks,Automatic}},GridLines->{Table[i,{i,Length[mean]}],Automatic},PlotRange->{{0,Length[mean]+1},plotRange},PlotStyle->{Automatic,Opacity[0],Opacity[0]},PlotTheme->OptionValue[PlotTheme]]]


(* ::Input::Initialization:: *)
plotRadianceVsRadiance[data_,OptionsPattern[{BinSize->{1.0},PlotRange->{{5,155},{5,155}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[data],OptionValue[BinSize],{"Log","PDF"},FrameLabel->{"Radiance ([radiance])","Radiance ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotResidualsVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-5,5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","K residual ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-7,7}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","Normalised K residual"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotResidualsHistogram[data_,OptionsPattern[{BinSize->{0.05},PlotRange->Automatic,PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{hist,plot,pr,mean,sdev,median,mdev},hist=Histogram[data,OptionValue[BinSize],"PDF",FrameLabel->{"K residual ([radiance])","Probability density"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]];
mean=Mean[data];sdev=StandardDeviation[data];median=Median[data];mdev=MeanDeviation[data];
pr=PlotRange/.Options[hist,PlotRange];
plot=Plot[{PDF[NormalDistribution[mean,sdev],x],PDF[LaplaceDistribution[median,mdev],x]},{x,pr[[1,1]],pr[[1,2]]},PlotRange->pr,PlotTheme->{"Monochrome","ThinLines"}];
Show[hist,plot]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsHistogram[data_,OptionsPattern[{BinSize->{0.1},PlotRange->Automatic,PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{hist,plot,pr,mean,sdev,median,mdev},hist=Histogram[data,OptionValue[BinSize],"PDF",FrameLabel->{"Normalised K residual","Probability density"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]];
mean=Mean[data];sdev=StandardDeviation[data];median=Median[data];mdev=MeanDeviation[data];
pr=PlotRange/.Options[hist,PlotRange];
plot=Plot[{PDF[NormalDistribution[mean,sdev],x],PDF[LaplaceDistribution[median,mdev],x]},{x,pr[[1,1]],pr[[1,2]]},PlotRange->pr,PlotTheme->{"Monochrome","ThinLines"}];
Show[hist,plot]];
