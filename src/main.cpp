#include "Processor.hpp"
#include <core/Simulation.hpp>
#include <boost/program_options.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/mpi/environment.hpp>

using namespace std;
using namespace sim;
namespace opts = boost::program_options;

#ifndef DIM
#define DIM 2
#endif

int main(int argc, char* argv[])
{
	/*
	 * Setup command line args
	 */
	opts::positional_options_description pcmd;
	pcmd.add("root_in",1);
	pcmd.add("root_out",1);
	pcmd.add("procs",1);
	pcmd.add("start",1);
	pcmd.add("stop",1);

	opts::options_description cmd1("Positonal Arguments");
	cmd1.add_options()("root_in",opts::value<string>(),"File name root for input .dat files.")
						 ("root_out",opts::value<string>(),"File name root for output .vtp files.")
						 ("procs",opts::value<int>(),"Number of processors.")
						 ("start",opts::value<int>(),"Timestep at which to start.")
						 ("stop",opts::value<int>(),"Timestep at which to stop.")
						 ;

	opts::options_description cmd2("Optional Commands");
	cmd2.add_options()("help","Print this message and exit.")
						 ("colours","Output particles from each colour separately")
						 ("iso","Generate an iso surface as well as the particle output (3D only).")
						 ("ip",opts::value<int>()->default_value(-1),"VALUE is the phase whose velocity, etc to interpolate onto the surface, if blank all phases are used.")
						 ("iv",opts::value<double>()->default_value(.25),"VALUE is the value of particle density at which to take the interface iso surface.")
						 ("pv",opts::value<double>()->default_value(.03),"VALUE is the attenuation coefficient ot use in generating the projections.")
						 ("nofail","If a file is not found, continue with next timestep rather than exiting.")
						 ("stride",opts::value<int>()->default_value(1),"Only process every VALUE-th file.")
						 ;

	opts::options_description cmd_opts("Command line arguments");
	cmd_opts.add(cmd1).add(cmd2);
	opts::variables_map cmd;
	opts::store(opts::command_line_parser(argc, argv).options(cmd_opts).positional(pcmd).style(opts::command_line_style::unix_style|opts::command_line_style::allow_long_disguise).run(), cmd);
	opts::notify(cmd);
	opts::arg = "VALUE";

	// print help if requested
	if(cmd.count("help"))
	{
		cout << "Usage: ubertrager ROOT PROCS START STOP [options]" << endl;
		cout << endl;
		cout << setprecision(3) << cmd_opts << endl;
		return 0;
	}

	/*
	 * Even though we don't use MPI, Simulation<DIM> objects use
	 * MPI routines in the constructor so we must initialize MPI.
	 */
	boost::mpi::environment mpi_env(argc,argv);

	// create the objects to read into
	size_t nprocs = cmd["procs"].as<int>();
	vector<Simulation<DIM>> cores(nprocs);

	cout << "nprocs = " << nprocs << endl;

	for(int t=cmd["start"].as<int>(); t<=cmd["stop"].as<int>(); t+=cmd["stride"].as<int>())
	{
		// aggregates and outputs the particles.
		Processor<DIM,2,2> processor;

		cout << "Loading timestep " << t << ": ";

		for(size_t i=0; i<nprocs; ++i)
		{
			cout << i << ", ";

			// construct file name
			stringstream fname;
			fname << cmd["root_in"].as<string>() << "." << t << "." << i << ".dat";

			// open file
			ifstream fin(fname.str());
			if(!fin.is_open())
			{
				stringstream msg;
				msg << "Error opening file " << fname.str() << "!";
				if(!cmd.count("nofail"))
				{
					throw runtime_error(msg.str());
				}
				else
				{
					cout << msg.str() << endl;
					continue;
				}
			}

			// deserialize
			boost::archive::binary_iarchive iarchive(fin);
			iarchive >> cores[i];

			// add the particles to the processor
			for(auto part : cores[i].fluidParticles())
				processor.addParticle(part);

			for(auto part : cores[i].wallParticles())
				processor.addParticle(part);
		}
		cout << endl;

		// generate vtk particles
		processor.generateVTKPointData(-1);

		stringstream fout_name;
		fout_name << cmd["root_out"].as<string>() << "." << t << ".vtp";
		cout << "Writing to " << fout_name.str() << endl;
		processor.write_vtp_file(fout_name.str());
	}

}
