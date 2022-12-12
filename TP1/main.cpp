#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/geometry/curves/GreedySegmentation.h>
#include "DGtal/io/Color.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

using namespace std;
using namespace DGtal;
using namespace Z2i;

#define M_PI 3.14159265358979323846

template<class T>
Curve getBoundary(T & object)
{
    KSpace kSpace; // Khalimsky space
    // We need to add a margine to prevent situations such that an object touch the bourder of the domain
    kSpace.init( object.domain().lowerBound() - Point(1,1), object.domain().upperBound() + Point(1,1), true);

    // 1) Call Surfaces::findABel() to find a cell which belongs to the border
    std::vector<Z2i::Point> boundaryPoints; // Boundary points are going to be stored here
    SCell aCell = Surfaces<Z2i::KSpace>::findABel(kSpace, object.pointSet(), 10000);

    // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object
    Curve boundaryCurve;
    SurfelAdjacency<2> SAdj( true );
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(boundaryPoints, kSpace, SAdj, object.pointSet(), aCell);

    // 3) Create a curve from a vector
    boundaryCurve.initFromVector(boundaryPoints);
    return boundaryCurve;
}

template<class T>
void sendToBoard( Board2D & board, T & p_Object, DGtal::Color p_Color) {
    board << CustomStyle( p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
}

int main(int argc, char** argv)
{
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet; // Digital set type
    //typedef Object<DT8_4, DigitalSet> ObjectType; // Digital object type
    typedef Object<DT4_8, DigitalSet> ObjectType; // Digital object type

    Board2D aBoard, aBoard2, aBoard3, aBoard3Line; // we will use these objects to save output
    Image image = PGMReader<Image>::importPGM (argv[1]); // you have to provide a correct path as a parameter
    
    // 1) Create a digital set of proper size
    DigitalSet set2d (image.domain());

    // 2) Use SetFromImage::append() to populate a digital set from the input image
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 0, 255);

    // 3) Create a digital object from the digital set    
    ObjectType set( dt4_8, set2d );
    std::vector< ObjectType > objects; // All conected components are going to be stored in it
    std::back_insert_iterator< std::vector< ObjectType > > inserter( objects ); // Iterator used to populated "objects".
    
    // 4) Use method writeComponents to obtain connected components
    set.writeComponents(inserter);
    std::cout << "Number of components : " << objects.size() << " (before elimination). \n"; // Right now size of "objects" is the number of conected components

    // 4,2) Elimination border grains
    objects.erase(std::remove_if(objects.begin(), objects.end(),
                                [image](ObjectType o) -> bool  
                                { 
                                    for(auto p : o.pointSet()) {
                                        if(p[0] == image.domain().lowerBound()[0] || 
                                           p[0] == image.domain().upperBound()[0] ||
                                           p[1] == image.domain().lowerBound()[1] || 
                                           p[1] == image.domain().upperBound()[1])
                                           { return true; }
                                    };
                                    return false;
                                }), objects.end());
    std::cout << "Number of components : " << objects.size() << " without border grains. \n"; // Right now size of "objects" is the number of conected components

    // 5 (Step 3) boundary
    Curve bondary = getBoundary(objects[0]);
    aBoard2 << bondary; 

    // 6 (Step 4)
    typedef Curve::PointsRange Range; // Container of digital points
    typedef Range::ConstIterator ConstIterator; // Iterator on the container
    typedef StandardDSS4Computer<ConstIterator> SegmentComputer; // StandardDSS4 computer
    typedef GreedySegmentation<SegmentComputer> Segmentation;

    Range range = bondary.getPointsRange(); // Construction of the computer
    // Segmentation
    SegmentComputer recognitionAlgorithm;
    Segmentation theSegmentation(range.begin(), range.end(), recognitionAlgorithm);

    for ( Segmentation::SegmentComputerIterator it = theSegmentation.begin(), itEnd = theSegmentation.end(); it != itEnd; ++it ) 
    {     
        // For ..
        aBoard3 << SetMode( "ArithmeticalDSS", "Points" ) << it->primitive(); // Cube line (optional)
        aBoard3 << SetMode( "ArithmeticalDSS", "BoundingBox" ) // Blue rectangle
                << CustomStyle( "ArithmeticalDSS/BoundingBox", new CustomPenColor( Color::Blue ) )
                << it->primitive();

        // For draw outEPS_Step4-2 image, with line
        aBoard3Line.drawLine(it->back()[0], it->back()[1], it->front()[0], it->front()[1]);
    } 

    float v1_x, v1_y, v2_x, v2_y, area = 0., perimeter = 0.;
    vector<float> sum_area_1, sum_area_2, sum_perimeter_1, sum_perimeter_2, sum_circularity_1, sum_circularity_2;
    // 7 (Step 4.2-5)
    for(auto o : objects) 
    { 
        Curve b = getBoundary(o); 
        Range r = b.getPointsRange();
        SegmentComputer recognitionAlgorithm;
        Segmentation theSegmentation(r.begin(), r.end(), recognitionAlgorithm);
        auto it0 = theSegmentation.begin();
        area = 0.;
        perimeter = sqrt(pow(it0->front()[0] - it0->back()[0] , 2) + pow(it0->front()[1] - it0->back()[1], 2));

        for ( Segmentation::SegmentComputerIterator it = ++it0, itEnd = theSegmentation.end(); it != itEnd; ++it ) 
        {     
            v1_x = it->back()[0]  - it0->back()[0];
            v1_y = it->back()[1]  - it0->back()[1];
            v2_x = it->front()[0] - it0->back()[0];
            v2_y = it->front()[1] - it0->back()[1];
            area += (v1_x * v2_y - v2_x * v1_y);
            perimeter += sqrt(pow(it->front()[0] - it->back()[0] , 2) + pow(it->front()[1] - it->back()[1], 2));
        } 
        area /= 2.;
        sum_area_1.emplace_back(o.size());
        sum_area_2.emplace_back(area);
        sum_perimeter_1.emplace_back(b.size());
        sum_perimeter_2.emplace_back(perimeter);
        sum_circularity_1.emplace_back((4. * M_PI * area) / (perimeter * perimeter));
        sum_circularity_2.emplace_back((4. * M_PI * o.size()) / (b.size() * b.size()));
    }

    std::cout << "--------------------------------------------\nAverage of the area of all grains of rice: \n";
    std::cout << " - with number of grid points : " << accumulate(sum_area_1.begin(), sum_area_1.end(), 0.0) / sum_area_1.size() << endl;
    std::cout << " - with area of the polygon : " << accumulate(sum_area_2.begin(), sum_area_2.end(), 0.0) / sum_area_2.size() << endl;
    std::cout << "--------------------------------------------\n";
    std::cout << "Average of the perimeter of all grains of rice: \n";
    std::cout << " - with the boundary  : " << accumulate(sum_perimeter_1.begin(), sum_perimeter_1.end(), 0.0) / sum_perimeter_1.size() << endl;
    std::cout << " - with the perimeter of the polygon  : " << accumulate(sum_perimeter_2.begin(), sum_perimeter_2.end(), 0.0) / sum_perimeter_2.size() << endl;
    std::cout << "--------------------------------------------\n";
    std::cout << "Average of the circularity of all grains of rice: \n";
    std::cout << " - with the cells  : " << accumulate(sum_circularity_1.begin(), sum_circularity_1.end(), 0.0) / sum_circularity_1.size() << endl;
    std::cout << " - with the polygon  : " << accumulate(sum_circularity_2.begin(), sum_circularity_2.end(), 0.0) / sum_circularity_2.size() << endl;

    // Save files
    sendToBoard(aBoard, objects[0], Color::Red); // Use this function to save digital objects to a file
    aBoard.saveEPS("outEPS_Step2.eps"); 
    aBoard2.saveEPS("outEPS_Step3.eps");
    aBoard3.saveEPS("outEPS_Step4.eps");
    aBoard3Line.saveEPS("outEPS_Step4-2.eps");
    return 0;
}
