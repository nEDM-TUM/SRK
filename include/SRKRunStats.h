#ifndef SRKRUNSTATS_H_
#define SRKRUNSTATS_H_

////////////////////////////////////////////////////////////////
/// SRKRunStats
///
/// The summarized results of a simulation run
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

struct SRKRunStats
{
	int numEvents = 0;  //Number of events
	double sXDetProb = 0.; //Cumulative probability of measuring up in the z direction

	double phiMean = 0.;  //Mean phi
	double phiError = 0.; //Uncertainty on Phi
	double phiStDev = 0.; //Standard deviation
	double phiKurtosis = 0.; //Kurtosis
	double phiKurtosisError = 0.; //Uncertainty on kurtosis
	double phiSkewness = 0.; //Skewness
	double phiSkewnessError = 0.; //Uncertainty on skewness
	double phiTsallisPower = 0.; //Tsallis power
	double phiTsallisPowerError = 0.; //Uncertainty on Tsallis power

	double thetaMean = 0.; //See phi definitions.
	double thetaError = 0.;
	double thetaStDev = 0.;
	double thetaKurtosis = 0.;
	double thetaKurtosisError = 0.;
	double thetaSkewness = 0.;
	double thetaSkewnessError = 0.;
	double thetaTsallisPower = 0.;
	double thetaTsallisPowerError = 0.;

};

#endif /* SRKRUNSTATS_H_ */
