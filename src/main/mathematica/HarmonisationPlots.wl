(* ::Package:: *)

(* ::Title:: *)
(*A Wolfram Language Package to Create Harmonisation Plots with Mathematica*)


(* ::Section::Closed:: *)
(*License*)


(* ::Text:: *)
(*Copyright (C) 2017 Ralf Quast*)
(**)
(*This code was developed for the EC project "Fidelity and Uncertainty in*)
(*Climate Data Records from Earth Observations (FIDUCEO)".*)
(*Grant Agreement: 638822*)
(*<Version> Reviewed and approved by <name, instituton>, <date>*)
(**)
(*This program is free software; you can redistribute it and/or modify it*)
(*under the terms of the GNU General Public License as published by the Free*)
(*Software Foundation; either version 3 of the License, or (at your option)*)
(*any later version.*)
(**)
(*This program is distributed in the hope that it will be useful, but WITHOUT*)
(*ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or*)
(*FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for*)
(*more details.*)
(**)
(*A copy of the GNU General Public License should have been supplied along*)
(*with this program; if not, see http://www.gnu.org/licenses/*)


(* ::Section:: *)
(*Plot functions*)


(* ::Input::Initialization:: *)
plotCorrelationMatrix[data_,OptionsPattern[{ColorFunction->"TemperatureMap",Frame->True,FrameLabel->{"Sensor","Sensor"},FrameTicks->None,Mesh->None,ParameterGrouping->0,PlotTheme->"Detailed"}]]:=Module[{labels,matrix,frameTicks,mesh,parameterGrouping},matrix=data[[1]];labels=data[[2]];
parameterGrouping=OptionValue[ParameterGrouping];frameTicks=If[parameterGrouping>0,Table[Table[{(parameterGrouping+1)/2+i-1,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {1,}"]->""]),45 Degree]},{i,1,Length[labels],parameterGrouping}],2],OptionValue[FrameTicks]];
mesh=If[parameterGrouping>0,{Length[matrix]/parameterGrouping-1,Length[matrix]/parameterGrouping-1},OptionValue[Mesh]];MatrixPlot[matrix,Mesh->mesh,ColorFunction->OptionValue[ColorFunction],Frame->OptionValue[Frame],FrameLabel->OptionValue[FrameLabel],FrameTicks->frameTicks,PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotCovarianceMatrix[data_,OptionsPattern[{ColorFunction->"TemperatureMap",Frame->True,FrameLabel->{"Sensor","Sensor"},FrameTicks->None,Mesh->None,ParameterGrouping->0,PlotTheme->"Detailed"}]]:=Module[{labels,matrix,frameTicks,mesh,parameterGrouping},matrix=data[[1]];labels=data[[2]];
parameterGrouping=OptionValue[ParameterGrouping];frameTicks=If[parameterGrouping>0,Table[Table[{(parameterGrouping+1)/2+i-1,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {1,}"]->""]),45 Degree]},{i,1,Length[labels],parameterGrouping}],2],OptionValue[FrameTicks]];
mesh=If[parameterGrouping>0,{Length[matrix]/parameterGrouping-1,Length[matrix]/parameterGrouping-1},OptionValue[Mesh]];MatrixPlot[matrix,Mesh->mesh,ColorFunction->OptionValue[ColorFunction],Frame->OptionValue[Frame],FrameLabel->OptionValue[FrameLabel],FrameTicks->frameTicks,PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotMeanResiduals[data_,OptionsPattern[{Exaggerate->100,PlotRange->Automatic,PlotTheme->{"Detailed","OpenMarkersThick","VibrantColor"}}]]:=Module[{exaggerate,plotRange,mean,stdev,range,labels,frameTicks},mean=data[[1]];exaggerate=OptionValue[Exaggerate];
plotRange=OptionValue[PlotRange];stdev=data[[2]];range=Max[Abs[{mean-exaggerate stdev,mean+exaggerate stdev}]];labels=data[[3]];frameTicks=Table[{i,Rotate[(StringReplace[StringJoin[labels[[i]]],RegularExpression[" {2,}"]->""]),90 Degree]},{i,Length[labels]}];ListPlot[{mean,mean-exaggerate stdev,mean+exaggerate stdev},Filling->{1->{2},1->{3}},FillingStyle->Thick,FrameLabel->{"Sensors","Mean K residual ([radiance])"},FrameTicks->{{Automatic,Automatic},{frameTicks,Automatic}},GridLines->{Table[i,{i,Length[mean]}],Automatic},PlotRange->{{0,Length[mean]+1},plotRange},PlotStyle->{Automatic,Opacity[0],Opacity[0]},PlotTheme->OptionValue[PlotTheme]]]


(* ::Input::Initialization:: *)
plotRadianceVsRadiance[data_,OptionsPattern[{BinSize->{1.0},PlotRange->{{5,155},{5,155}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[data],OptionValue[BinSize],{"Log","PDF"},FrameLabel->{"Radiance ([radiance])","Radiance ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotRadianceVsRadianceRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{{5,155},{5,155}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[data],OptionValue[PointNumber]],FrameLabel->{"Radiance ([radiance])","Radiance ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotRadianceVsRadianceRandomLegended[data_,names_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{{5,155},{5,155}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{randomData,pointNumbers},RandomSeed[OptionValue[Seed]];
pointNumbers=Table[Round[OptionValue[PointNumber]Length[data[[i,1]]]/Sum[Length[data[[j,1]]],{j,Length[data]}]],{i,Length[data]}];
randomData=Table[RandomChoice[Transpose[{data[[i,1]],data[[i,2]]}],pointNumbers[[i]]],{i,Length[data]}];
randomData=Table[Legended[randomData[[i]],names[[i]]],{i,Length[data]}];ListPlot[randomData,FrameLabel->{"Radiance ([radiance])","Radiance ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotRadianceVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{1.0}},PlotRange->{All,{-5,5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","Radiance ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotRadianceVsTimeRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{5,155}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],FrameLabel->{"Time (calendar year)","Radiance ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotRadianceUncertaintyVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-5,5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","Radiance uncertainty ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotRadianceUncertaintyVsTimeRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,Automatic},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],FrameLabel->{"Time (calendar year)","Radiance uncertainty ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotRadianceDifferenceVsRadiance[data_,OptionsPattern[{BinSize->{1.0,0.1},PlotRange->{{5,155},{-5,-5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[data],OptionValue[BinSize],{"Log","PDF"},FrameLabel->{"Radiance ([radiance])","Radiance difference ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotRadianceDifferenceVsRadianceRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{{5,155},{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[data],OptionValue[PointNumber]],FrameLabel->{"Radiance ([radiance])","Radiance difference ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotRadianceDifferenceVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-5,5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","Radiance difference ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotRadianceDifferenceVsTimeRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],FrameLabel->{"Time (calendar year)","Radiance difference ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotResidualsVsRadiance[data_,OptionsPattern[{BinSize->{1.0,0.1},PlotRange->{{5,155},{-5,-5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[data],OptionValue[BinSize],{"Log","PDF"},FrameLabel->{"Radiance ([radiance])","K residual ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotResidualsVsRadianceRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{{5,155},{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[data],OptionValue[PointNumber]],FrameLabel->{"Radiance ([radiance])","K residual ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotResidualsVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-5,5}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","K residual ([radiance])"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotResidualsVsTimeRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],FrameLabel->{"Time (calendar year)","K residual ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotResidualsVsTimeRandomLegended[data_,names_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{randomData,pointNumbers},RandomSeed[OptionValue[Seed]];
pointNumbers=Table[Round[OptionValue[PointNumber]Length[data[[i,1]]]/Sum[Length[data[[j,1]]],{j,Length[data]}]],{i,Length[data]}];
randomData=Table[RandomChoice[Transpose[{data[[i,1]]/86400/365.25+1970,data[[i,2]]}],pointNumbers[[i]]],{i,Length[data]}];
randomData=Table[Legended[randomData[[i]],names[[i]]],{i,Length[data]}];ListPlot[randomData,FrameLabel->{"Time (calendar year)","K residual ([radiance])"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsVsTime[data_,OptionsPattern[{BinSizes->{{1/4},{0.1}},PlotRange->{All,{-7,7}},PlotTheme->"Detailed",Title->None}]]:=DensityHistogram[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[BinSizes],{"Log","PDF"},FrameLabel->{"Time (calendar year)","Normalised K residual"},ColorFunction->"TemperatureMap",PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotRangePadding->None,PlotTheme->OptionValue[PlotTheme]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsVsTimeRandom[data_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{},RandomSeed[OptionValue[Seed]];ListPlot[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],FrameLabel->{"Time (calendar year)","Normalised K residual"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsVsTimeRandomLegended[data_,names_,OptionsPattern[{Seed->27181,PointNumber->10000,PlotRange->{All,{-5,5}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Module[{randomData,pointNumbers},RandomSeed[OptionValue[Seed]];
pointNumbers=Table[Round[OptionValue[PointNumber]Length[data[[i,1]]]/Sum[Length[data[[j,1]]],{j,Length[data]}]],{i,Length[data]}];
randomData=Table[RandomChoice[Transpose[{data[[i,1]]/86400/365.25+1970,data[[i,2]]}],pointNumbers[[i]]],{i,Length[data]}];
randomData=Table[Legended[randomData[[i]],names[[i]]],{i,Length[data]}];ListPlot[randomData,FrameLabel->{"Time (calendar year)","Normalised K residual"},PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotTheme->OptionValue[PlotTheme]]];


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


(* ::Input::Initialization:: *)
computeTrend[data_,OptionsPattern[{ConfidenceLevel->0.99}]]:=Module[{lm},lm=LinearModelFit[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],{1,t},t];lm[{"BestFit","MeanPredictionBands"},ConfidenceLevel->OptionValue[ConfidenceLevel]]];


(* ::Input::Initialization:: *)
computeTrendRandom[data_,OptionsPattern[{ConfidenceLevel->0.99,Seed->27181,PointNumber->100000}]]:=Module[{lm},RandomSeed[OptionValue[Seed]];lm=LinearModelFit[RandomChoice[Transpose[{data[[1]]/86400/365.25+1970,data[[2]]}],OptionValue[PointNumber]],{1,t},t];lm[{"BestFit","MeanPredictionBands"},ConfidenceLevel->OptionValue[ConfidenceLevel]]];


(* ::Input::Initialization:: *)
plotNormalizedResidualsTrend[trend_,tmin_,tmax_,OptionsPattern[{GridLines->Automatic,PlotRange->{All,{-0.05,0.05}},PlotTheme->{"Detailed","VibrantColor"},Title->None}]]:=Plot[Evaluate[Flatten[trend]],{t,tmin,tmax},Axes->False,Frame->True,FrameLabel->{"Time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (calendar year)","Mean normalised K residual"},GridLines->OptionValue[GridLines],PlotLabel->OptionValue[Title],PlotRange->OptionValue[PlotRange],PlotStyle->{Automatic,{Black,Dashed,Thin},{Black,Dashed,Thin}},PlotTheme->OptionValue[PlotTheme],PlotLabels->{Placed[trend[[1]],{Scaled[1],Before}],None,None},PlotLegends->None];
