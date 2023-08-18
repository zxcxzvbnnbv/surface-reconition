#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCurvatures.h>
#include <vtkSTLWriter.h>
#include <vtkPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkThreshold.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataArray.h>
#include <vtkAutoInit.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkQuadricDecimation.h>
#include <vtkIdList.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <vector>
#include <set>
#include <ctime>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);
using namespace std;

enum surfaceType {
	Peak,
	Ridge,
	Saddle_ridge,
	notExist,
	Flat,
	Minimal,
	Pit,
	Valley,
	Saddle_valley 
};
const double thresholdK =  4e-4;
const double thresholdH = 1e-4;

void AddLabel(vtkPolyData* polyData);
void AddPointLabel(vtkPolyData* polyData);
void findNeighbors( vtkPolyData* polyData, vtkIdType seedID, vtkIdList* neighbors);
void findPointNeighbors(vtkPolyData* polyData, vtkIdType seedID, vtkIdList* pointNeighbors);
bool fusion(vtkPolyData* polyData, vtkIntArray* pointLabels, int threshold, int* iteration);
void RegionGrow(vtkPolyData* polyData, double thresholdK, double thresholdH, vtkIntArray* pointLabels,  vtkIntArray* triangleLabels);
void regionGrow( vtkPolyData* polyData, double thresholdK, double thresholdH, vtkIntArray* labels, vtkIntArray* surfType);
surfaceType surfTypedef(double gaussCurvature, double meanCurvature, double thresholdK, double thresholdH);
int main()
{
	
	// 读取STL文件
	vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
	stlReader->SetFileName("E:\\application\\opencascade\\blade2.stl");
	stlReader->Update();

	vtkPolyData* polyData = stlReader->GetOutput();

    // Create a label array to store region labels for each triangle
	vtkSmartPointer<vtkIntArray> pointLabels = vtkSmartPointer<vtkIntArray>::New();
	pointLabels->SetName("PointRegionLabels");
	pointLabels->SetNumberOfComponents(1);
	pointLabels->SetNumberOfTuples(polyData->GetNumberOfPoints());
	pointLabels->FillComponent(0, -1); // Initialize with default value

	vtkSmartPointer<vtkIntArray> triangleLabels = vtkSmartPointer<vtkIntArray>::New();
	triangleLabels->SetName("TriangleRegionLabels");
	triangleLabels->SetNumberOfComponents(1);
	triangleLabels->SetNumberOfTuples(polyData->GetNumberOfCells());
	triangleLabels->FillComponent(0, -1); // Initialize with default value

	vtkSmartPointer<vtkIntArray> surfType = vtkSmartPointer<vtkIntArray>::New();
	surfType->SetName("surfaceTypeLabels");
	surfType->SetNumberOfComponents(1);
	surfType->SetNumberOfTuples(polyData->GetNumberOfCells());
	surfType->FillComponent(0, -1); // Initialize with default value
	
	//三角片曲率类型判定方法一
	//AddLabel(polyData);
	//regionGrow(polyData, thresholdK, thresholdH, triangleLabels, surfType);

	//三角片曲率类型判定方法二
	AddPointLabel(polyData);
    RegionGrow(polyData, thresholdK, thresholdH, pointLabels, triangleLabels);

	/*
		for (vtkIdType i = 0; i < pointLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = pointLabels->GetValue(i);
		std::cout << "第" << i << "个点的标签为：" << currentValue << std::endl;
	}

	int maxValue = triangleLabels->GetValue(0);  // 假设数组至少有一个元素
		for (vtkIdType i = 0; i < triangleLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = triangleLabels->GetValue(i);
		if (currentValue > maxValue)
		{
			maxValue = currentValue;
		}
	}

	int LabelNumber[2000] = { 0 };
	for (int i = 0; i < maxValue; ++i)
	{
		int numberOfLabel = 0;
		for (vtkIdType j = 0; j < polyData->GetNumberOfCells(); ++j)
		{
			if (triangleLabels->GetValue(j) == i) {
				numberOfLabel++;
			}
		}
		LabelNumber[i] = numberOfLabel;
	}

	for (int i = 0; i < maxValue; ++i)
	{
		std::cout << "第" << i << "个数据块有" << LabelNumber[i] << "个三角片" << std::endl;
	}
	*/


	int maxValue = triangleLabels->GetValue(0);  // 假设数组至少有一个元素
	for (vtkIdType i = 0; i < triangleLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = triangleLabels->GetValue(i);
		if (currentValue > maxValue)
		{
			maxValue = currentValue;
		}
	}
	
	


	// 获取三角形的数量
	vtkCellArray* triangles = polyData->GetPolys();
	vtkIdType ntriangles = triangles->GetNumberOfCells();

	// 创建颜色数组，用于为每个三角形设置颜色
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	// 设置三角形颜色（这里示意性地设置了随机颜色）
	for (vtkIdType i = 0; i < ntriangles; ++i)
	{
		double color[3] = { 255,255,255 };
		if (triangleLabels->GetValue(i) != -1)
		{
			double rgb = 255 * (maxValue - triangleLabels->GetValue(i)) / maxValue;
			color[0] = rgb;
			color[1] = std::abs(rgb - 100);
			color[2] = std::abs(rgb * 2 - 255);
		}
		/*		  
			if (triangleLabels->GetValue(i) == 228)
			{
				color[1] = 0;
				color[2] = 0;
			}
		    double rgb = 255 * (maxValue - triangleLabels->GetValue(i)) / maxValue;
			color[0] = rgb;
			color[1] = std::abs(rgb - 100);
			color[2] = std::abs(rgb * 2 - 255); 
		*/

		colors->InsertTuple(i, color);
	}

	polyData->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkCellArray> nTriangles = vtkSmartPointer<vtkCellArray>::New();
	for (vtkIdType cellId = 0; cellId < polyData->GetNumberOfCells(); ++cellId)
	{
		if (triangleLabels->GetValue(cellId) == 0)
		{
			   vtkCell* cell = polyData->GetCell(cellId);
			   nTriangles->InsertNextCell(cell);
		}	
	}

	vtkSmartPointer<vtkPolyData> npolyData = vtkSmartPointer<vtkPolyData>::New();
	npolyData->SetPolys(nTriangles);

	//新的mapper和actor
	vtkSmartPointer<vtkPolyDataMapper> partMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	partMapper->SetInputData(npolyData);

	vtkSmartPointer<vtkActor> partActor = vtkSmartPointer<vtkActor>::New();
	partActor->SetMapper(partMapper);
	
	//mapper和actor
	vtkSmartPointer<vtkPolyDataMapper> fullMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	fullMapper->SetInputData(polyData);

	vtkSmartPointer<vtkActor> fullActor = vtkSmartPointer<vtkActor>::New();
	fullActor->SetMapper(fullMapper);

	// 创建渲染器和渲染窗口
	double fullView[4] = { 0,0,1,1 };
	double partView[4] = { 0.55,0,1,1 };
	vtkSmartPointer<vtkRenderer> fullRenderer = vtkSmartPointer<vtkRenderer>::New();
	fullRenderer->SetViewport(fullView);
    fullRenderer->AddActor(fullActor);

	vtkSmartPointer<vtkRenderer> partRenderer = vtkSmartPointer<vtkRenderer>::New();
	partRenderer->SetViewport(partView);
	partRenderer->AddActor(partActor);


	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(fullRenderer);
	//renderWindow->AddRenderer(partRenderer);
	renderWindow->SetSize(960, 320);
    renderWindow->Render();

	// 创建交互器
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();


	return 0;
}




//计算模型的高斯和平均曲率作为标签
void AddLabel(vtkPolyData* polyData)
{

	// 获取曲率数据
	vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
	curvaturesFilter->SetInputData(polyData);
	curvaturesFilter->SetCurvatureTypeToMean();
	curvaturesFilter->Update();
	vtkDataArray* pmeanCurvatureData = curvaturesFilter->GetOutput()->GetPointData()->GetArray("Mean_Curvature");
	pmeanCurvatureData->SetName("MeanCurvature");
	

	vtkSmartPointer<vtkCurvatures> curvaturesFilter2 = vtkSmartPointer<vtkCurvatures>::New();
	curvaturesFilter2->SetInputData(polyData);
	curvaturesFilter2->SetCurvatureTypeToGaussian();
	curvaturesFilter2->Update();
	vtkDataArray* pgaussCurvatureData = curvaturesFilter2->GetOutput()->GetPointData()->GetArray("Gauss_Curvature");
	pgaussCurvatureData->SetName("GaussianCurvature");
	

	// 遍历每个三角片，获取高斯曲率和平均曲率
	vtkIdType numCells = polyData->GetNumberOfCells();
	vtkSmartPointer<vtkDoubleArray> gaussCurvatureData = vtkSmartPointer<vtkDoubleArray>::New();
	gaussCurvatureData->SetName("GaussianCurvature");
	gaussCurvatureData->SetNumberOfComponents(1);
	gaussCurvatureData->SetNumberOfTuples(numCells);
	
	vtkSmartPointer<vtkDoubleArray> meanCurvatureData = vtkSmartPointer<vtkDoubleArray>::New();
	meanCurvatureData->SetName("MeanCurvature");
	meanCurvatureData->SetNumberOfComponents(1);
	meanCurvatureData->SetNumberOfTuples(numCells);

	
	for (vtkIdType cellId = 0; cellId < numCells; ++cellId)
	{
		double gaussCurvature, meanCurvature;
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(polyData->GetCell(cellId));
		if (triangle)
		{
			gaussCurvature = pgaussCurvatureData->GetTuple1(triangle->GetPointId(0));
			meanCurvature = pmeanCurvatureData->GetTuple1(triangle->GetPointId(0));
			gaussCurvatureData->SetValue(cellId, gaussCurvature);
			meanCurvatureData->SetValue(cellId, meanCurvature);
		}
	}
	// Attach curvature arrays to the polydata's cell data
	polyData->GetCellData()->AddArray(gaussCurvatureData);
	polyData->GetCellData()->AddArray(meanCurvatureData);

}


void AddPointLabel(vtkPolyData* polyData)
{

	// 获取曲率数据
	vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
	curvaturesFilter->SetInputData(polyData);
	curvaturesFilter->SetCurvatureTypeToMean();
	curvaturesFilter->Update();
	vtkDataArray* pmeanCurvatureData = curvaturesFilter->GetOutput()->GetPointData()->GetArray("Mean_Curvature");
	pmeanCurvatureData->SetName("MeanCurvature");


	vtkSmartPointer<vtkCurvatures> curvaturesFilter2 = vtkSmartPointer<vtkCurvatures>::New();
	curvaturesFilter2->SetInputData(polyData);
	curvaturesFilter2->SetCurvatureTypeToGaussian();
	curvaturesFilter2->Update();
	vtkDataArray* pgaussCurvatureData = curvaturesFilter2->GetOutput()->GetPointData()->GetArray("Gauss_Curvature");
	pgaussCurvatureData->SetName("GaussianCurvature");

	polyData->GetPointData()->AddArray(pmeanCurvatureData);
	polyData->GetPointData()->AddArray(pgaussCurvatureData);
}

//寻找当前三角片相邻的三角片：计算复杂度太高，待优化
void findNeighbors(vtkPolyData* polyData, vtkIdType seedID, vtkIdList* neighbors)
{
	vtkCellArray* triangles = polyData->GetPolys();
	vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
	vtkIdType sp1, sp2, sp3, p1, p2, p3;
	double spoint[3][3];

	// The current triangle is adjacent to the seed triangle
	vtkTriangle* seedTriangle = dynamic_cast<vtkTriangle*>(polyData->GetCell(seedID));

	if (seedTriangle->GetCellType() == VTK_TRIANGLE)
	{
		sp1 = seedTriangle->GetPointId(0);
		sp2 = seedTriangle->GetPointId(1);
		sp3 = seedTriangle->GetPointId(2);

		points->GetPoint(sp1, spoint[0]);
		points->GetPoint(sp2, spoint[1]);
		points->GetPoint(sp3, spoint[2]);
	}

	/*
	std::cout << "Point 1: " << spoint[0][0] << " " << spoint[0][1] << " " << spoint[0][2] << std::endl;
	std::cout << "Point 2: " << spoint[1][0] << " " << spoint[1][1] << " " << spoint[1][2] << std::endl;
	std::cout << "Point 3: " << spoint[2][0] << " " << spoint[2][1] << " " << spoint[2][2] << std::endl;
	std::cout << "------------------------" << std::endl;
	*/


	for (vtkIdType triangleId = 0; triangleId < polyData->GetNumberOfCells(); ++triangleId)
	{
		if (triangleId == seedID) continue;
		//通过查找是否有两个相同顶点判断是否为相邻三角形
		int sharedPointsCount = 0;
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(polyData->GetCell(triangleId));
		double point[3][3];
		p1 = triangle->GetPointId(0);
		p2 = triangle->GetPointId(1);
		p3 = triangle->GetPointId(2);

		points->GetPoint(p1, point[0]);
		points->GetPoint(p2, point[1]);
		points->GetPoint(p3, point[2]);

		for (vtkIdType i = 0; i < 3; ++i)
		{
			for (vtkIdType j = 0; j < 3; ++j)
			{
				if ((spoint[i][0]==point[j][0])&& (spoint[i][1] == point[j][1]) && (spoint[i][2] == point[j][2]))
				{
					sharedPointsCount++;
					break;
				}
			}
		}
		if (sharedPointsCount >= 2)
		{
			// The current triangle is adjacent to the seed triangle
			neighbors->InsertNextId(triangleId);
		}
	}

}

//判断三角片类型
surfaceType surfTypedef(double gaussCurvature, double meanCurvature, double thresholdK, double thresholdH)
{
	if (gaussCurvature > thresholdK && meanCurvature > thresholdH)   return surfaceType::Peak;
	if (thresholdK > gaussCurvature && gaussCurvature > -thresholdK && meanCurvature > thresholdH)   return surfaceType::Ridge;
	if (gaussCurvature < -thresholdK && meanCurvature > thresholdH)   return surfaceType::Saddle_ridge;
	if (gaussCurvature > thresholdK && thresholdH > meanCurvature > -thresholdH)   return surfaceType::notExist;
	if (thresholdK > gaussCurvature && gaussCurvature > -thresholdK &&  thresholdH > meanCurvature > -thresholdH)   return surfaceType::Flat;
	if (gaussCurvature < -thresholdK && thresholdH > meanCurvature > -thresholdH)   return surfaceType::Minimal;
	if (gaussCurvature > thresholdK && meanCurvature < -thresholdH)   return surfaceType::Pit;
	if (thresholdK > gaussCurvature && gaussCurvature > -thresholdK && meanCurvature < -thresholdH)   return surfaceType::Valley;
	if (gaussCurvature < -thresholdK && meanCurvature < -thresholdH)   return surfaceType::Saddle_valley;
}

//区域生长方法一：计算复杂度太高，待优化
void regionGrow(vtkPolyData* polyData, double thresholdK, double thresholdH, vtkIntArray* labels, vtkIntArray* surfType)
{
	vtkDataArray* meanCurvatureData = polyData->GetCellData()->GetArray("MeanCurvature");
	vtkDataArray* gaussCurvatureData = polyData->GetCellData()->GetArray("GaussianCurvature");
	std::queue<vtkIdType> regionQueue;
	vtkIdType cellnum = polyData->GetNumberOfCells();
	vtkIdType seedId = 0;
	surfaceType blockType, faceType;

	double gaussCur, meanCur, gaussCurvature, meanCurvature;

	int blockId = 0;
	for (vtkIdType i = 0; i < cellnum; ++i)
	{  
        
		//确定块初始生长点的位置和类型
		seedId = i;
		if (labels->GetValue(seedId) != -1) continue;
		gaussCurvature = gaussCurvatureData->GetTuple1(seedId);
		meanCurvature = meanCurvatureData->GetTuple1(seedId);
		blockType = surfTypedef(gaussCurvature, meanCurvature, thresholdK, thresholdH);
		surfType->SetVariantValue(seedId, blockType);
		labels->SetValue(seedId, blockId);
		regionQueue.push(seedId);

		//块区域生长
		while (!regionQueue.empty())
		{
			vtkIdList* neighbors = vtkIdList::New();
			vtkIdType triangelId = regionQueue.front();
			regionQueue.pop();
			findNeighbors(polyData, triangelId, neighbors);
			for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); ++j)
			{
				vtkIdType id = neighbors->GetId(j);
				// 处理 id
				if (labels->GetValue(id) == -1)
				{
					gaussCur = gaussCurvatureData->GetTuple1(id);
					meanCur = meanCurvatureData->GetTuple1(id);
					faceType = surfTypedef(gaussCur, meanCur, thresholdK, thresholdH);
					if (faceType == blockType)
					{
						regionQueue.push(id);
						surfType->SetValue(id, faceType);
						labels->SetValue(id, blockId);
					}
				}
			}
			neighbors->Delete();
		}
	 //更新块
			  ++blockId;
	}
}

/*
--------------------------------------------------------------------------
寻找点的相邻点：
先找到包含该点的所有三角片，再把三角片上的所有点放入集合
--------------------------------------------------------------------------
*/
void findPointNeighbors(vtkPolyData* polyData, vtkIdType seedID, vtkIdList* pointNeighbors)
{
	vtkNew<vtkIdList>  cellIds;
	polyData->GetPointCells(seedID, cellIds);
	std::set<vtkIdType> neighbours;
	vtkIdType cellNum = cellIds->GetNumberOfIds();
	for (vtkIdType i = 0; i < cellNum; ++i)
	{
		auto cellId = cellIds->GetId(i);
		vtkNew<vtkIdList> cellPointIds;
		polyData->GetCellPoints(cellId, cellPointIds);
		for (vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); ++j)
		{
			neighbours.insert(cellPointIds->GetId(j));
		}
	}
	for(std::set<vtkIdType>::iterator it = neighbours.begin(); it != neighbours.end(); it++)
	{
		pointNeighbors->InsertNextId(*it);
	}
	
}

/*
--------------------------------------------------------------------------
区域生长方法二：
先根据点的曲率对点进行划分，再根据点所属的数据块判断三角
片所属的数据块
--------------------------------------------------------------------------
*/
void RegionGrow(vtkPolyData* polyData, double thresholdK, double thresholdH, vtkIntArray* pointLabels, vtkIntArray* triangleLabels)
{
	vtkDataArray* meanCurvatureData = polyData->GetPointData()->GetArray("MeanCurvature");
	vtkDataArray* gaussCurvatureData = polyData->GetPointData()->GetArray("GaussianCurvature");
	std::queue<vtkIdType> regionQueue;
	vtkIdType pointNum = polyData->GetNumberOfPoints();
	vtkIdType seedPointId = 0;
	surfaceType blockType, pointType;

	double gaussCur, meanCur, gaussCurvature, meanCurvature;

	int blockId = 0;
	for (vtkIdType i = 0; i < pointNum; ++i)
	{
         //确定种子点的位置和类型
		seedPointId = i;
		if (pointLabels->GetValue(seedPointId) != -1) continue;
		gaussCurvature = gaussCurvatureData->GetTuple1(seedPointId);
		meanCurvature = meanCurvatureData->GetTuple1(seedPointId);
		blockType = surfTypedef(gaussCurvature, meanCurvature, thresholdK, thresholdH);
		pointLabels->SetValue(seedPointId, blockId);
		regionQueue.push(seedPointId);

		//块区域生长
		while (!regionQueue.empty())
		{
			vtkIdList* pointNeighbors = vtkIdList::New();
			vtkIdType pointId = regionQueue.front();
			regionQueue.pop();
			findPointNeighbors(polyData, pointId, pointNeighbors);
			vtkIdType pointNum = pointNeighbors->GetNumberOfIds();
			for (vtkIdType j = 0; j < pointNum; ++j)
			{
				vtkIdType id = pointNeighbors->GetId(j);
				// 处理 id
				if (pointLabels->GetValue(id) == -1)
				{
					gaussCur = gaussCurvatureData->GetTuple1(id);
					meanCur = meanCurvatureData->GetTuple1(id);
					pointType = surfTypedef(gaussCur, meanCur, thresholdK, thresholdH);
					if (pointType == blockType)
					{
						regionQueue.push(id);
						pointLabels->SetValue(id, blockId);
					}
				}
			}
			pointNeighbors->Delete();
		}
		//更新块
		++blockId;
	}

	int iteration = 0;
    fusion(polyData, pointLabels, 5, &iteration);
	
	/*
	int maxValue = pointLabels->GetValue(0);  // 假设数组至少有一个元素
	for (vtkIdType i = 0; i < pointLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = pointLabels->GetValue(i);
		if (currentValue > maxValue)
		{
			maxValue = currentValue;
		}
	}
	int LabelNumber[2000] = { 0 };
	for (int i = 0; i < maxValue; ++i)
	{
		int numberOfLabel = 0;
		for (vtkIdType j = 0; j < polyData->GetNumberOfPoints(); ++j)
		{
			if (pointLabels->GetValue(j) == i) {
				numberOfLabel++;
			}
		}
		LabelNumber[i] = numberOfLabel;
	}

	for (int i = 0; i < maxValue; ++i)
	{
		std::cout << "第" << i << "个数据块有" << LabelNumber[i] << "个三角片" << std::endl;
	}
	*/
	
	



	//根据点判断面所属的块
	vtkIdType numCell = polyData->GetNumberOfCells();
	for (vtkIdType cellId = 0; cellId < numCell; ++cellId)
	{
		auto seed = time(0);
		srand(seed);
		vtkNew<vtkIdList> cellPointIds;
		polyData->GetCellPoints(cellId, cellPointIds);
		vtkIdType id[4];
		for (vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); ++j)
		{
			id[j] = cellPointIds->GetId(j);
		}

		if (pointLabels->GetValue(id[0]) == pointLabels->GetValue(id[1])
			&& pointLabels->GetValue(id[1]) == pointLabels->GetValue(id[2]))
		{
			triangleLabels->SetValue(cellId, pointLabels->GetValue(id[0]));
			continue;
		}

		if (pointLabels->GetValue(id[0]) == pointLabels->GetValue(id[1]))
		{
			triangleLabels->SetValue(cellId, pointLabels->GetValue(id[0]));
		}

		if (pointLabels->GetValue(id[0]) == pointLabels->GetValue(id[2]))
		{
			triangleLabels->SetValue(cellId, pointLabels->GetValue(id[0]));
		}

		if (pointLabels->GetValue(id[1]) == pointLabels->GetValue(id[2]))
		{
			triangleLabels->SetValue(cellId, pointLabels->GetValue(id[1]));
		}
		else {
			triangleLabels->SetValue(cellId, pointLabels->GetValue(id[rand()%3]));
		}

	}

}

/*
--------------------------------------------------------------------------
数据块融合：
根据阈值，将三角片数量较少的块合并到临近的较大的数据块中
若没有低于阈值的块，返回0，反之返回1；
--------------------------------------------------------------------------
*/
bool fusion(vtkPolyData* polyData, vtkIntArray* pointLabels, int threshold, int* iteration)
{
	if (*iteration > 5)  return false;
	int maxValue = pointLabels->GetValue(0);  // 假设数组至少有一个元素
	for (vtkIdType i = 0; i < pointLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = pointLabels->GetValue(i);
		if (currentValue > maxValue)
		{
			maxValue = currentValue;
		}
	}

	vector<int> LabelNumber;
	vector<vector<int>> LabelBlock;
	
	for (int i = 0; i < maxValue; ++i)
	{
		int numberOfLabel = 0;
		vector<int> LabelIndex;
		for (vtkIdType j = 0; j < polyData->GetNumberOfPoints(); ++j)
		{
			if (pointLabels->GetValue(j) == i) {
				numberOfLabel++;
				LabelIndex.push_back(j);
			}
		}
		if (numberOfLabel != 0)
		{
		   LabelBlock.push_back(LabelIndex);
		   LabelNumber.push_back(numberOfLabel);
		}


	}

	//若没有低于阈值的块，返回0
	int minValue = *min_element(LabelNumber.begin(), LabelNumber.end());
	if (minValue > threshold)  return false;

	//查找数量低于阈值的块，进行合并
	for (int i = 0; i < LabelNumber.size(); ++i)
	{
		if (LabelNumber[i] < threshold)
		{
			//先找块之外的所有临近点
			set<vtkIdType> neighbours;
			for (int j = 0; j < LabelNumber[i]; ++j)
			{
				vtkNew<vtkIdList>  cellIds;
				polyData->GetPointCells(LabelBlock[i][j], cellIds);
				vtkIdType cellNum = cellIds->GetNumberOfIds();
				for (vtkIdType cellIndex = 0; cellIndex < cellNum; ++cellIndex)
				{
					auto cellId = cellIds->GetId(cellIndex);
					vtkNew<vtkIdList> cellPointIds;
					polyData->GetCellPoints(cellId, cellPointIds);
					for (vtkIdType pointIndex = 0; pointIndex < cellPointIds->GetNumberOfIds(); ++pointIndex)
					{
						neighbours.insert(cellPointIds->GetId(pointIndex));
					}
				}
			}

			//找到数量最大的临近点所属的块，
			//如果周围的临近点所属块数量小于原本块则不变
			//否则，将原来的点合并到新的块里

			int maxLabelNum = 0;
			int maxLabel = 0; 
			int num = 0;
			for (set<vtkIdType>::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
			{
				num = LabelNumber[pointLabels->GetValue(*it)];
				if (num > maxLabelNum)
				{
					maxLabelNum = num;
					maxLabel = pointLabels->GetValue(*it);
				}
			}
			for (int m = 0; m < LabelNumber[i]; ++m)
			{
				pointLabels->SetValue(LabelBlock[i][m], maxLabel);
				LabelBlock[maxLabel].push_back(LabelBlock[i][m]);
			}
			sort(LabelBlock[maxLabel].begin(), LabelBlock[maxLabel].end());
			LabelNumber[maxLabel] += LabelNumber[i];
			LabelNumber[i] = 0;
			LabelBlock[i].clear();
		}
	}
	++(*iteration);
	return true;
}

/*
bool fusion(vtkPolyData* polyData, vtkIntArray* pointLabels, int threshold, int* iteration)
{
	if (*iteration > 5)  return false;
	int maxValue = pointLabels->GetValue(0);  // 假设数组至少有一个元素
	for (vtkIdType i = 0; i < pointLabels->GetNumberOfTuples(); ++i)
	{
		int currentValue = pointLabels->GetValue(i);
		if (currentValue > maxValue)
		{
			maxValue = currentValue;
		}
	}

	vector<int> LabelNumber;
	vector<vector<int>> LabelBlock;

	for (int i = 0; i < maxValue; ++i)
	{
		int numberOfLabel = 0;
		vector<int> LabelIndex;
		for (vtkIdType j = 0; j < polyData->GetNumberOfPoints(); ++j)
		{
			if (pointLabels->GetValue(j) == i) {
				numberOfLabel++;
				LabelIndex.push_back(j);
			}
		}
		if (numberOfLabel != 0)
		{
		   LabelBlock.push_back(LabelIndex);
		   LabelNumber.push_back(numberOfLabel);
		}


	}

	//若没有低于阈值的块，返回0
	int minValue = *min_element(LabelNumber.begin(), LabelNumber.end());
	if (minValue > threshold)  return false;

	//查找数量低于阈值的块，进行合并
	for (int i = 0; i < LabelNumber.size(); ++i)
	{
		if (LabelNumber[i] < threshold)
		{
			//先找块之外的所有临近点
			set<vtkIdType> neighbours;
			for (int j = 0; j < LabelNumber[i]; ++j)
			{
				vtkNew<vtkIdList>  cellIds;
				polyData->GetPointCells(LabelBlock[i][j], cellIds);
				vtkIdType cellNum = cellIds->GetNumberOfIds();
				for (vtkIdType cellIndex = 0; cellIndex < cellNum; ++cellIndex)
				{
					auto cellId = cellIds->GetId(cellIndex);
					vtkNew<vtkIdList> cellPointIds;
					polyData->GetCellPoints(cellId, cellPointIds);
					for (vtkIdType pointIndex = 0; pointIndex < cellPointIds->GetNumberOfIds(); ++pointIndex)
					{
						neighbours.insert(cellPointIds->GetId(pointIndex));
					}
				}
			}

			//找到数量最大的临近点所属的块，
			//如果周围的临近点所属块数量小于原本块则不变
			//否则，将原来的点合并到新的块里

			int maxLabelNum = 0;
			int maxLabel = 0;
			int maxIndex = 0;
			int num = 0;
			int actuIndex = 0;
			for (set<vtkIdType>::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
			{
				for (int index = 0; index < LabelBlock.size(); ++index)
				{
					vector<int> labellist = LabelBlock[index];
					auto t = find(labellist.begin(), labellist.end(), *it);
					if (t != labellist.end())
					{
						num = LabelNumber[index];
						actuIndex = index;
						break;
					}
				}
				if (num > maxLabelNum)
				{
					maxLabelNum = num;
					maxLabel = pointLabels->GetValue(*it);
					maxIndex = actuIndex;
				}
			}
			for (int m = 0; m < LabelNumber[i]; ++m)
			{
				pointLabels->SetValue(LabelBlock[i][m], maxLabel);
				LabelBlock[maxIndex].push_back(LabelBlock[i][m]);
			}
			sort(LabelBlock[maxIndex].begin(), LabelBlock[maxIndex].end());
			LabelNumber[maxIndex] += LabelNumber[i];
			LabelNumber[i] = 0;
			LabelBlock[i].clear();
		}
	}
	++(*iteration);
	return true;
}
*/