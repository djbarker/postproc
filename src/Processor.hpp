#ifndef PROCESSOR_HPP_
#define PROCESSOR_HPP_

#include <vector>
#include <Particle.hpp>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetAttributes.h>
#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkImageExtractComponents.h>

using namespace std;
using namespace sim;

template<size_t D, size_t T, size_t C>
class Processor
{
public:
	Processor(){};
	virtual ~Processor(){};

	void addParticle(sim::Particle<D,T,C> p);

	void load_dat_file(std::string fname, bool root);
	void write_vtp_file(std::string fname);
	void initialise_sublists();
	void remesh(int vel_fluid = -1);
	void generateVTKPointData(int fluid);
	void generateIsoSurface(int fluid, double iso_val);
	void generateMesh();

protected:

	void write_vtkPolyData(std::string);
	void write_vtkImageData(std::string);

	vector<Particle<D,T,C>> particles;

	vtkSmartPointer<vtkMarchingCubes> surface;
	vtkSmartPointer<vtkImageData> voxels;
	vtkSmartPointer<vtkImageData> fluid_voxels;		// store separately as we don't really want this in the output, just to isosurface it.
	vtkSmartPointer<vtkImageExtractComponents> ext_col;
	vtkSmartPointer<vtkProbeFilter> probe;
	vtkSmartPointer<vtkDoubleArray> rel_density;
	vtkSmartPointer<vtkDoubleArray> pressure;
	vtkSmartPointer<vtkDoubleArray> velocity;
	vtkSmartPointer<vtkDoubleArray> part_id;
	vtkSmartPointer<vtkDataSet> vtp_output; // data to actually write to file
};

template<size_t D, size_t T, size_t C>
void Processor<D,T,C>::addParticle(Particle<D,T,C> p)
{
    particles.push_back(p);
}

template<size_t D, size_t T, size_t C>
void Processor<D,T,C>::generateVTKPointData(int fluid)
{
    const int tout = 0; // which timestep to output

	vtkSmartPointer<vtkPolyData> dataSet = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIdList> IdList = vtkSmartPointer<vtkIdList>::New();

	newPts->SetDataTypeToDouble();
	int num = 0;

	for(auto part : particles)
	{
		if ((fluid<0) || (part.fluid==fluid))
		{
			if(part.type!=GhostP)
			{
			if(D==3) newPts->InsertNextPoint(discard_dims(part.pos[tout][0]),
									discard_dims(part.pos[tout][1]),
									discard_dims(part.pos[tout][2]) );
			else     newPts->InsertNextPoint(discard_dims(part.pos[tout][0]),
									discard_dims(part.pos[tout][1]),
									0.0 );
			}
			else     newPts->InsertNextPoint(discard_dims(part.pos[tout][0]),
												discard_dims(part.pos[tout][1]),
												0.1 );
			IdList->InsertNextId(num);
			num++;
		}
	}

	dataSet->Allocate(num);
	dataSet->SetPoints(newPts);
	dataSet->Modified();
	dataSet->Update();

	dataSet->Allocate(num);
	dataSet->InsertNextCell(VTK_POLY_VERTEX, IdList);

	vtkSmartPointer<vtkDoubleArray> PressureArray = vtkSmartPointer<vtkDoubleArray>::New();
	PressureArray->SetName("Pressure");
	PressureArray->SetNumberOfComponents(1);
	PressureArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkIntArray> ColourArray = vtkSmartPointer<vtkIntArray>::New();
	ColourArray->SetName("Fluid");
	ColourArray->SetNumberOfComponents(1);
	ColourArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkDoubleArray> DensityArray = vtkSmartPointer<vtkDoubleArray>::New();
	DensityArray->SetName("Density");
	DensityArray->SetNumberOfComponents(1);
	DensityArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkDoubleArray> SigmaArray = vtkSmartPointer<vtkDoubleArray>::New();
	SigmaArray->SetName("Sigma");
	SigmaArray->SetNumberOfComponents(1);
	SigmaArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkDoubleArray> VelocityArray = vtkSmartPointer<vtkDoubleArray>::New();
	VelocityArray->SetName("Velocity");
	VelocityArray->SetNumberOfComponents(3);
	VelocityArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkDoubleArray> ForceArray = vtkSmartPointer<vtkDoubleArray>::New();
	ForceArray->SetName("Acceleration");
	ForceArray->SetNumberOfComponents(3);
	ForceArray->SetNumberOfTuples(num);

	vtkSmartPointer<vtkIntArray> IdArray = vtkSmartPointer<vtkIntArray>::New();
	IdArray->SetName("Id");
	IdArray->SetNumberOfComponents(1);
	IdArray->SetNumberOfTuples(num);

	num=0;
	for(auto& part : particles)
	{
		if ((fluid<0) || ((int)part.fluid==fluid))
		{
			PressureArray->InsertValue(num, discard_dims(part.pressure));
			ColourArray->InsertValue(num, part.fluid);
			IdArray->InsertValue(num, part.id);
			DensityArray->InsertValue(num, discard_dims(part.density[tout]));

			//double vx = (part.v[0][0]+part.v[1][0])*0.5;
			//double vy = (part.v[0][1]+part.v[1][1])*0.5;
			//double vz = (part.v[0][2]+part.v[1][2])*0.5;
			//VelocityArray->InsertTuple3(num,vx,vy,vz);
			VelocityArray->InsertTuple3(num, discard_dims(part.vel[tout][0]),
											 discard_dims(part.vel[tout][1]),
											 (D==3?discard_dims(part.vel[tout][2]):0.0) );

			ForceArray->InsertTuple3(num, discard_dims(part.acc[0]),
										  discard_dims(part.acc[1]),
										  (D==3?discard_dims(part.acc[2]):0.0) );

            SigmaArray->InsertValue(num,discard_dims(part.sigma));
			//if(Particle::sa_boundaries)	GammaArray->InsertValue(num,part.gamma[0]);

			num++;
		}
	}

	if (fluid>=0)
	{
		cout << "Total particles of fluid " << fluid << ": " << num << "\t(" << 100*(double)num/(double)particles.size() << "%)" << endl;
	}
	else if (fluid == -1)
	{
		cout << "Total particles: " << num << endl;
	}

	dataSet->GetPointData()->AddArray(PressureArray);
	dataSet->GetPointData()->AddArray(ColourArray);
	dataSet->GetPointData()->SetActiveAttribute("Colour",vtkDataSetAttributes::SCALARS);
	dataSet->GetPointData()->AddArray(IdArray);
	dataSet->GetPointData()->AddArray(DensityArray);
	dataSet->GetPointData()->AddArray(VelocityArray);
	dataSet->GetPointData()->AddArray(ForceArray);
	dataSet->GetPointData()->AddArray(SigmaArray);
	//if(Particle::sa_boundaries)	dataSet->GetPointData()->AddArray(GammaArray);

	dataSet->Modified();
	dataSet->Update();

	vtp_output = dataSet;
}

template<size_t D, size_t T, size_t C>
void Processor<D,T,C>::write_vtp_file(string fname)
{
	if(vtp_output->IsA("vtkPolyData"))
	{
		write_vtkPolyData(fname);
	}
	else
	{
		write_vtkImageData(fname);
	}
}


template<size_t D, size_t T, size_t C>
void Processor<D,T,C>::write_vtkPolyData(string fname)
{
	vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
	vtkSmartPointer<vtkXMLPolyDataWriter> writer= vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetDataModeToBinary();
	writer->SetFileName(fname.c_str());
	writer->SetInput(vtp_output);
	writer->SetCompressor(compressor);

	if (!writer->Write())
	{
		throw runtime_error("Error writing polydata to vtk file!");
	}
}

template<size_t D, size_t T, size_t C>
void Processor<D,T,C>::write_vtkImageData(string fname)
{
	vtkSmartPointer<vtkZLibDataCompressor> compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();
	vtkSmartPointer<vtkXMLImageDataWriter> writer= vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writer->SetDataModeToBinary();
	writer->SetFileName(fname.c_str());
	writer->SetInput(vtp_output);
	writer->SetCompressor(compressor);

	if (!writer->Write())
	{
		throw runtime_error("Error writing imagedata to vtk file!");
	}
}

#endif /* PROCESSOR_HPP_ */
