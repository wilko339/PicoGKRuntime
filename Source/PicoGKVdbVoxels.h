//
// SPDX-License-Identifier: Apache-2.0
//
// PicoGK ("peacock") is a compact software kernel for computational geometry,
// specifically for use in Computational Engineering Models (CEM).
//
// For more information, please visit https://picogk.org
//
// PicoGK is developed and maintained by LEAP 71 - © 2023 by LEAP 71
// https://leap71.com
//
// Computational Engineering will profoundly change our physical world in the
// years ahead. Thank you for being part of the journey.
//
// We have developed this library to be used widely, for both commercial and
// non-commercial projects alike. Therefore, have released it under a permissive
// open-source license.
//
// The foundation of PicoGK is a thin layer on top of the powerful open-source
// OpenVDB project, which in turn uses many other Free and Open Source Software
// libraries. We are grateful to be able to stand on the shoulders of giants.
//
// LEAP 71 licenses this file to you under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with the
// License. You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, THE SOFTWARE IS
// PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.
//
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef PICOGKVDBVOXELS_H_
#define PICOGKVDBVOXELS_H_

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/RayIntersector.h>

#include "PicoGKMesh.h"

using namespace openvdb;

#define PICOGK_VOXEL_DEFAULTBACKGROUND 3.0f

namespace PicoGK
{

class Voxels
{

public:
    typedef std::shared_ptr<Voxels> Ptr;
    
    Voxels(float fBackground = PICOGK_VOXEL_DEFAULTBACKGROUND)
    {
        m_roFloatGrid = FloatGrid::create(fBackground);
        m_roFloatGrid->setGridClass(GRID_LEVEL_SET);
    };
    
    Voxels( float fBackground,
            FloatGrid::Ptr roGrid)
    {
        m_roFloatGrid = roGrid;
    };
    
    Voxels(const Voxels& oSource)
    {
        m_roFloatGrid   = deepCopyTypedGrid<FloatGrid>(oSource.m_roFloatGrid);
    };

    ~Voxels()
    {
    }
    
    bool bIsEqual(const Voxels& oCompare) const
    {
        CoordBBox oBBoxThis = m_roFloatGrid->evalActiveVoxelBoundingBox();
        CoordBBox oBBoxComp = oCompare.m_roFloatGrid->evalActiveVoxelBoundingBox();
        
        int32_t iMinX = std::min(   oBBoxThis.min().x(),
                                    oBBoxComp.min().x());
        
        int32_t iMinY = std::min(   oBBoxThis.min().y(),
                                    oBBoxComp.min().y());
        
        int32_t iMinZ = std::min(   oBBoxThis.min().z(),
                                    oBBoxComp.min().z());
        
        int32_t iMaxX = std::max(   oBBoxThis.max().x(),
                                    oBBoxComp.max().x());
        
        int32_t iMaxY = std::max(   oBBoxThis.max().y(),
                                    oBBoxComp.max().y());
        
        int32_t iMaxZ = std::max(   oBBoxThis.min().z(),
                                    oBBoxComp.max().z());
        
        auto oThis = m_roFloatGrid->getConstAccessor();
        auto oComp = oCompare.m_roFloatGrid->getConstAccessor();
        
        for(int32_t x = iMinX; x <= iMaxX; x++)
        for(int32_t y = iMinY; y <= iMaxY; y++)
        for(int32_t z = iMinZ; z <= iMaxZ; z++)
        {
            openvdb::Coord xyz(x,y,z);
            
            bool bThisInside = (oThis.getValue(xyz) <= 0.0f);
            bool bCompInside = (oComp.getValue(xyz) <= 0.0f);
            
            if (bThisInside != bCompInside)
                return false;
        }
        
        return true;
     }

    void BoolAdd(const Voxels& oOther)
    {
        FloatGrid::Ptr roOperand = deepCopyTypedGrid<FloatGrid>(oOther.m_roFloatGrid);
        openvdb::tools::csgUnion(*m_roFloatGrid, *roOperand);
    }

    void BoolSubtract(const Voxels& oOther)
    {
        FloatGrid::Ptr roOperand = deepCopyTypedGrid<FloatGrid>(oOther.m_roFloatGrid);
        openvdb::tools::csgDifference(*m_roFloatGrid, *roOperand);
    }

    void BoolIntersect(const Voxels& oOther)
    {
        FloatGrid::Ptr roOperand = deepCopyTypedGrid<FloatGrid>(oOther.m_roFloatGrid);
        openvdb::tools::csgIntersection(*m_roFloatGrid, *roOperand);
    }
    
    void Offset(float fSize, VoxelSize oVoxelSize)
    {
        openvdb::tools::LevelSetFilter<openvdb::FloatGrid> oFilter(*m_roFloatGrid);
        
        float fSizeVx = -oVoxelSize.fToVoxels(fSize); // openvdb treats offsets as inwards
        
        // The following command doesn't seem to be necessary, but verify
        //oFilter.resize(std::abs(fSizeVx) + fBackground());
        oFilter.offset(fSizeVx);
    }
    
    void DoubleOffset(  float fSize1,
                        float fSize2,
                        VoxelSize oVoxelSize)
    {
        openvdb::tools::LevelSetFilter<openvdb::FloatGrid> oFilter(*m_roFloatGrid);
        
        float fSize1Vx = -oVoxelSize.fToVoxels(fSize1); // openvdb treats offsets as inwards
        float fSize2Vx = -oVoxelSize.fToVoxels(fSize2); // openvdb treats offsets as inwards
        
        oFilter.offset(fSize1Vx);
        oFilter.offset(fSize2Vx);
    }
    
    void TripleOffset(  float fSize,
                        VoxelSize oVoxelSize)
    {
        openvdb::tools::LevelSetFilter<openvdb::FloatGrid> oFilter(*m_roFloatGrid);
        
        float fSizeVx = -oVoxelSize.fToVoxels(fSize); // openvdb treats offsets as inwards
        
        // offset inwards first
        oFilter.offset(-fSizeVx);
        
        // offset twice the size outwards next
        oFilter.offset(fSizeVx * 2);
        
        // offset inwards again. Now we are back where we started
        // but have lost a lot of detail = smooth
        oFilter.offset(-fSizeVx);
    }

    void RenderMesh(	const Mesh& oMesh,
    					VoxelSize oVoxelSize)
    {
        // We have to convert the mesh to voxel coords before
        // we transfer it to openvdb for rendering
        // We should use the openvdb transformations in the
        // future and inform the openvdb float grid of our
        // voxel size. This is, however, a bit of a larger
        // project and will be done at a later time.
        
        Mesh oMeshInVoxelCoord;
        for (int32_t n=0; n<oMesh.nTriangleCount(); n++)
        {
            Vector3 vecA(0,0,0);
            Vector3 vecB(0,0,0);
            Vector3 vecC(0,0,0);
            oMesh.GetTriangle(n, &vecA, &vecB, &vecC);
            
            vecA = oVoxelSize.vecToVoxels(vecA);
            vecB = oVoxelSize.vecToVoxels(vecB);
            vecC = oVoxelSize.vecToVoxels(vecC);
            
            oMeshInVoxelCoord.nAddTriangle(vecA, vecB, vecC);
        }
        
        FloatGrid::Ptr roVoxelized = roFloatGridFromMesh(   oMeshInVoxelCoord,
                                                            1.0f,
                                                            fBackground());
        
        openvdb::tools::csgUnion(*m_roFloatGrid, *roVoxelized);
    }
    
    void RenderLattice( const Lattice& oLattice,
                        float fVoxelSizeMM)
    {
        auto oAccess = m_roFloatGrid->getAccessor();
        
        for (auto roSphere : oLattice.oSpheres())
        {
            DoRenderLattice(&oAccess, fBackground(), *roSphere, fVoxelSizeMM);
        }
        
        for (auto roBeam : oLattice.oBeams())
        {
            DoRenderLattice(&oAccess, fBackground(), *roBeam, fVoxelSizeMM);
        }
    }
    
    void RenderImplicit(    const BBox3& oBBox,
                            PKPFnfSdf pfn,
                            VoxelSize oVoxelSize)
    {
        auto oAccess = m_roFloatGrid->getAccessor();
        
        Coord xyzMin = oVoxelSize.xyzToVoxels(oBBox.vecMin);
        Coord xyzMax = oVoxelSize.xyzToVoxels(oBBox.vecMax);
        
        // Increase the bounding box by the voxel distance of the background value
        // so we don't cut off the narrow band
        int32_t iAdd = (int32_t) (m_roFloatGrid->background() + 0.5f);
        
        for(int32_t x = xyzMin.X - iAdd; x <= xyzMax.X + iAdd; x++)
        for(int32_t y = xyzMin.Y - iAdd; y <= xyzMax.Y + iAdd; y++)
        for(int32_t z = xyzMin.Z - iAdd; z <= xyzMax.Z + iAdd; z++)
        {
            Vector3 vecSample = oVoxelSize.vecToMM(Coord(x,y,z));
            openvdb::Coord xyz(x,y,z);
            
            float fValue = std::min(    oVoxelSize.fToVoxels((*pfn)(&vecSample)),
                                        oAccess.getValue(xyz));
            
            if (fValue < fBackground())
            {
                // inside or near inside
                
                if (fValue > -fBackground())
                {
                    oAccess.setValue(xyz, fValue);
                }
                else
                {
                    oAccess.setValue(xyz, -fBackground());
                }
            }
        }
    }
    
    void IntersectImplicit( PKPFnfSdf pfn,
                            VoxelSize oVoxelSize)
    {
        Voxels oVox(fBackground());
        
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        
        BBox3 oBBoxMM;
        oBBoxMM.vecMin.X = oVoxelSize.fToMM(oBBox.min().x());
        oBBoxMM.vecMin.Y = oVoxelSize.fToMM(oBBox.min().y());
        oBBoxMM.vecMin.Z = oVoxelSize.fToMM(oBBox.min().z());
        
        oBBoxMM.vecMax.X = oVoxelSize.fToMM(oBBox.max().x());
        oBBoxMM.vecMax.Y = oVoxelSize.fToMM(oBBox.max().y());
        oBBoxMM.vecMax.Z = oVoxelSize.fToMM(oBBox.max().z());
        
        oVox.RenderImplicit(oBBoxMM, pfn, oVoxelSize);
        
        // Swap out the grids, so we keep using the "nice"
        // implict grid, and use our grid just as the mask
        // Not sure if this is really necessary, but the
        // implicit grid is "perfect"
        m_roFloatGrid.swap(oVox.m_roFloatGrid);
        
        BoolIntersect(oVox);
    }

    Mesh::Ptr roAsMesh(float fVoxelSizeMM, bool bTriangulateMeshes, float fMeshAdaptivity) const
    {
        Mesh::Ptr roMesh = std::make_shared<Mesh>();
        
    	std::vector< openvdb::Vec3s > oPoints;
    	std::vector< openvdb::Vec3I > oTriangles;
    	std::vector< openvdb::Vec4I > oQuads;
		
        openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*m_roFloatGrid,
                                                         oPoints,
                                                         oTriangles,
                                                         oQuads,
                                                         0.0f,
                                                         fMeshAdaptivity,
                                                         true);

        if (bTriangulateMeshes)
        {
            for (const openvdb::Vec4I oQuad : oQuads)
            {
                openvdb::Vec3I o1(oQuad[0], oQuad[1], oQuad[2]);
                oTriangles.push_back(o1);

                openvdb::Vec3I o2(oQuad[2], oQuad[3], oQuad[0]);
                oTriangles.push_back(o2);
            }
            oQuads.clear();
        }
	    
        for (const openvdb::Vec3s& v : oPoints)
        {
            Vector3 vec(v.x(), v.y(), v.z());
            vec *= fVoxelSizeMM;
            roMesh->nAddVertex(vec);
        }

        for (const openvdb::Vec3I oTri : oTriangles)
        {
            roMesh->nAddTriangle(   PicoGK::Triangle(   oTri[2],
                                                        oTri[1],
                                                        oTri[0]));
        }

        for (const openvdb::Vec4I oQuad : oQuads)
        {
            roMesh->nAddQuad(PicoGK::Quad(
                oQuad[3],
                oQuad[2],
                oQuad[1],
                oQuad[0]));
        }

        return roMesh;
    }
    
    void ProjectZSliceDn(   float fZStart,
                            float fZEnd,
                            VoxelSize oVoxelSize)
    {
        assert(fZStart > fZEnd);
        
        int32_t iZStart = oVoxelSize.iToVoxels(fZStart);
        int32_t iZEnd   = oVoxelSize.iToVoxels(fZEnd);
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        
        auto oAccess = m_roFloatGrid->getAccessor();
        
        for(int32_t x = oBBox.min().x(); x <= oBBox.max().x(); x++)
        for(int32_t y = oBBox.min().y(); y <= oBBox.max().y(); y++)
        for(int32_t z = iZStart; z > iZEnd; z--)
        {
            
            openvdb::Coord xyz(x,y,z);
            openvdb::Coord xyzUnder(x,y,z-1);
            
            float fValue = std::min(    oAccess.getValue(xyzUnder),
                                        oAccess.getValue(xyz));
            
            if (fValue < fBackground())
            {
                // inside or near inside
                
                if (fValue > -fBackground())
                {
                    oAccess.setValue(xyzUnder, fValue);
                }
                else
                {
                    oAccess.setValue(xyzUnder, -fBackground());
                }
            }
        }
    }
    
    void ProjectZSliceUp(   float fZStart,
                            float fZEnd,
                            VoxelSize oVoxelSize)
    {
        assert(fZStart < fZEnd);
        
        int32_t iZStart = oVoxelSize.iToVoxels(fZStart);
        int32_t iZEnd   = oVoxelSize.iToVoxels(fZEnd);
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        
        auto oAccess = m_roFloatGrid->getAccessor();
        
        for(int32_t x = oBBox.min().x(); x <= oBBox.max().x(); x++)
        for(int32_t y = oBBox.min().y(); y <= oBBox.max().y(); y++)
        for(int32_t z = iZStart; z <= iZEnd; z++)
        {
            openvdb::Coord xyz(x,y,z);
            openvdb::Coord xyzOver(x,y,z+1);
            
            float fValue = std::min(    oAccess.getValue(xyzOver),
                                        oAccess.getValue(xyz));
            
            if (fValue < fBackground())
            {
                // inside or near inside
                
                if (fValue > -fBackground())
                {
                    oAccess.setValue(xyzOver, fValue);
                }
                else
                {
                    oAccess.setValue(xyzOver, -fBackground());
                }
            }
        }
    }

    void ProjectZSlice( float fZStart,
                        float fZEnd,
                        VoxelSize oVoxelSize)
    {
        if (fZStart > fZEnd)
            ProjectZSliceDn(fZStart, fZEnd, oVoxelSize);
        else
            ProjectZSliceUp(fZStart, fZEnd, oVoxelSize);
    }
    
    void CalculateProperties(   float* pfVolume,
                                BBox3* poBBox,
                                VoxelSize oVoxelSize)
    {
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        
        auto oAccess = m_roFloatGrid->getConstAccessor();
        
        int nCount = 0;
        
        BBox3 oResult;
        
        for (int32_t x=oBBox.min().x(); x<=oBBox.max().x(); x++)
        for (int32_t y=oBBox.min().y(); y<=oBBox.max().y(); y++)
        for (int32_t z=oBBox.min().z(); z<=oBBox.max().z(); z++)
        {
            if (oAccess.getValue(openvdb::Coord(x,y,z)) <= 0.0f)
            {
                // Voxel is set
                nCount++;
                oResult.Include(oVoxelSize.vecToMM(Coord(x,y,z)));
            }
        }
        
        float fVolume = nCount;
        fVolume *= oVoxelSize;
        fVolume *= oVoxelSize;
        fVolume *= oVoxelSize; // cubic!
        
        *pfVolume   = fVolume;
        *poBBox     = oResult;
    }
    
    inline void GetSurfaceNormal(   Vector3 vecPt,
                                    VoxelSize oVoxelSize,
                                    Vector3* pvecNormal)
    {
        math::GradStencil oStencil(*m_roFloatGrid);
        Coord xyz = oVoxelSize.xyzToVoxels(vecPt);
    
        oStencil.moveTo(    openvdb::Coord( xyz.X,
                                            xyz.Y,
                                            xyz.Z));
        
        auto vecGradient = oStencil.gradient();
        vecGradient.normalize();
        
        pvecNormal->X = vecGradient.x();
        pvecNormal->Y = vecGradient.y();
        pvecNormal->Z = vecGradient.z();
    }
    
    inline bool bFindClosestPointOnSurface( Vector3 vecSearch,
                                            VoxelSize oVoxelSize,
                                            Vector3* pvecSurfacePoint)
    {
        Coord xyzSearch = oVoxelSize.xyzToVoxels(vecSearch);
        
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        oBBox.expand(openvdb::Coord(xyzSearch.X, xyzSearch.Y, xyzSearch.Z));
        // Add the search point to the BBox and expand a little to avoid
        // stopping too soon
        
        auto oAccess    = m_roFloatGrid->getConstAccessor();
        
        int32_t iMaxSearchRadius = (int32_t) std::ceil(sqrt(
                    oBBox.extents().x() * oBBox.extents().x() +
                    oBBox.extents().y() * oBBox.extents().y() +
                    oBBox.extents().z() * oBBox.extents().z()));
        
        bool bStartInside = oAccess.getValue(openvdb::Coord(xyzSearch.X, xyzSearch.Y, xyzSearch.Z)) <= 0.0f;
        
        for (int32_t iR=0; iR<iMaxSearchRadius; iR++)
        {
            bool bOutsideActiveBounds = true;
            Coord xyzSurfacePoint;
            
            if (bBresenhamSphereHitTest(    bStartInside,
                                            xyzSearch,
                                            iR,
                                            oAccess,
                                            oBBox,
                                            &xyzSurfacePoint,
                                            &bOutsideActiveBounds))
            {
                *pvecSurfacePoint = oVoxelSize.vecToMM(xyzSurfacePoint);
                return true;
            }
            
            if (bOutsideActiveBounds)
                return false;
        }
        
        return false;
    }
    
    inline bool bRayCastToSurface(  const Vector3& vecSearch,
                                    const Vector3& vecDirection,
                                    VoxelSize oVoxelSize,
                                    Vector3* pvecSurfacePoint)
    {
        tools::LevelSetRayIntersector oIntersector(*m_roFloatGrid);
        math::Ray<Real> oRay(   Vec3f(  oVoxelSize.fToVoxels(vecSearch.X),
                                        oVoxelSize.fToVoxels(vecSearch.Y),
                                        oVoxelSize.fToVoxels(vecSearch.Z)),
                                Vec3f(  vecDirection.X,
                                        vecDirection.Y,
                                        vecDirection.Z));
        
        math::Vec3<Real> xyz;
        math::Vec3<Real> vecNormal;
        
        if (oIntersector.intersectsIS(oRay, xyz))
        {
            *pvecSurfacePoint = oVoxelSize.vecToMM( Coord(  xyz.x(),
                                                            xyz.y(),
                                                            xyz.z()));
            return true;
        }

        // If the ray reaches its maximum range without hitting a filled voxel,
        // return false to indicate that there was no intersection.
            
        return false;
    }
    
    void GetVoxelDimensions(    int32_t* pnXSize,
                                int32_t* pnYSize,
                                int32_t* pnZSize) const
    {
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        *pnXSize = oBBox.extents().x();
        *pnYSize = oBBox.extents().y();
        *pnZSize = oBBox.extents().z();
    }
    
    void GetSlice( int32_t nZSlice,
                   float* pfBuffer)
    {
        CoordBBox oBBox = m_roFloatGrid->evalActiveVoxelBoundingBox();
        openvdb::Coord xyz(0, 0, nZSlice + oBBox.min().z());
        
        auto oAccess = m_roFloatGrid->getConstAccessor();
        
        int32_t n=0;
        for (xyz.y()=oBBox.min().y(); xyz.y()<=oBBox.max().y(); xyz.y()++)
        for (xyz.x()=oBBox.min().x(); xyz.x()<=oBBox.max().x(); xyz.x()++)
        {
            pfBuffer[n] = oAccess.getValue(xyz);
            n++;
        }
    }

    FloatGrid::Ptr roVdbGrid() const 	{return m_roFloatGrid;}
    
    inline float fBackground() const    {return m_roFloatGrid->background();}

    Coord 	        m_xyzSize;
    FloatGrid::Ptr	m_roFloatGrid;
    
    
protected:
    
    template<class TAccessor, class TLatticeBeam>
    static void DoRenderLattice(    TAccessor* poAccess,
                                    float fBackground,
                                    const TLatticeBeam& oLattice,
                                    VoxelSize oVoxelSize)
    {
        Vector3 vecMin(oLattice.vecMin());
        Vector3 vecMax(oLattice.vecMax());
            
        Coord xyzMin = oVoxelSize.xyzToVoxels(vecMin);
        Coord xyzMax = oVoxelSize.xyzToVoxels(vecMax);
        
        // Increase the bounding box by the voxel distance of the background value
        // so we don't cut off the narrow band
        int32_t iAdd = (int32_t) (fBackground + 0.5f);
        
        for(int32_t x = xyzMin.X - iAdd; x <= xyzMax.X + iAdd; x++)
        for(int32_t y = xyzMin.Y - iAdd; y <= xyzMax.Y + iAdd; y++)
        for(int32_t z = xyzMin.Z - iAdd; z <= xyzMax.Z + iAdd; z++)
        {
            openvdb::Coord xyz(x,y,z);
            
            Vector3 vecSample = oVoxelSize.vecToMM(Coord(x,y,z));
            
            // Boolean add to existing value, if one exists
            float fValue = std::min(    oVoxelSize.fToVoxels(oLattice.fSdValue(vecSample)),
                                        poAccess->getValue(xyz));
                    
            if (fValue < fBackground)
                poAccess->setValue(xyz, fValue);
        }
    }
    
    static FloatGrid::Ptr roFloatGridFromMesh(  const Mesh& oMesh,
                                                float fVoxelSizeMM,
                                                float fBackground)
    {
        std::vector< openvdb::Vec3s >   oVertices;
        std::vector< openvdb::Vec3I >   oTriangles;
        std::vector< openvdb::Vec4I >   oQuads;
        
        for (int32_t n=0; n<oMesh.nVertexCount(); n++)
        {
            Vector3 v(0.0f,0.0f,0.0f);
            oMesh.GetVertex(n, (Vector3*)&v);
            oVertices.push_back(openvdb::Vec3s( v.X,
                                                v.Y,
                                                v.Z));
        }
        
        // Purely quad mesh
        if (oMesh.nTriangleCount() == 0 && oMesh.nQuadCount() > 0)
        {
            for (int32_t n = 0; n < oMesh.nQuadCount(); n++)
            {
                Quad q;
                oMesh.GetQuad(n, &q);
                oQuads.push_back(openvdb::Vec4I(    (unsigned int)q.A,
                                                    (unsigned int)q.B,
                                                    (unsigned int)q.C,
                                                    (unsigned int)q.D));
            }

            openvdb::math::Transform::Ptr roTransform
                = openvdb::math::Transform::createLinearTransform(fVoxelSizeMM);

            return openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*roTransform,
                oVertices,
                oQuads,
                fBackground);
        }

        // Mixed quad and tri mesh (happens with Rhino)
        if (oMesh.nTriangleCount() > 0 && oMesh.nQuadCount() > 0)
        {
            for (int32_t n = 0; n < oMesh.nQuadCount(); n++)
            {
                Quad q;
                oMesh.GetQuad(n, &q);
                oQuads.push_back(openvdb::Vec4I((unsigned int)q.A,
                    (unsigned int)q.B,
                    (unsigned int)q.C,
                    (unsigned int)q.D));
            }

            for (int32_t n = 0; n < oMesh.nTriangleCount(); n++)
            {
                Triangle t;
                oMesh.GetTriangle(n, &t);
                oTriangles.push_back(openvdb::Vec3I((unsigned int)t.A,
                    (unsigned int)t.B,
                    (unsigned int)t.C));
            }

            openvdb::math::Transform::Ptr roTransform
                = openvdb::math::Transform::createLinearTransform(fVoxelSizeMM);

            return openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*roTransform,
                oVertices,
                oTriangles,
                oQuads,
                fBackground);
        }

        // The only remaining case should be a purely triangular mesh
        else
        {
            for (int32_t n=0; n<oMesh.nTriangleCount(); n++)
            {
                Triangle t;
                oMesh.GetTriangle(n, &t);
                oTriangles.push_back(openvdb::Vec3I(    (unsigned int) t.A,
                                                        (unsigned int) t.B,
                                                        (unsigned int) t.C));
            }

            openvdb::math::Transform::Ptr roTransform
                = openvdb::math::Transform::createLinearTransform(fVoxelSizeMM);

            return openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(  *roTransform,
                                                                        oVertices,
                                                                        oTriangles,
                                                                        fBackground);
        }
    }
    
    template <class TAccessor>
    inline bool bBresenhamSphereHitTest(    bool bReferenceInside,
                                            Coord xyzCenter,
                                            int iRadius,
                                            const TAccessor& oAccess,
                                            const openvdb::CoordBBox& oBB,
                                            Coord* pxyzPoint,
                                            bool* pbOutsideActiveBounds)
    {
        *pbOutsideActiveBounds = true;
        
        openvdb::Coord xyz;
        
        for (xyz.z() = xyzCenter.Z - iRadius; xyz.z() <= xyzCenter.Z + iRadius; xyz.z()++)
        {
            int r2 = iRadius * iRadius;
            for (xyz.y() = xyzCenter.Y - iRadius; xyz.y() <= xyzCenter.Y + iRadius; xyz.y()++)
            {
                for (xyz.x() = xyzCenter.X - iRadius; xyz.x() <= xyzCenter.X + iRadius; xyz.x()++)
                {
                    if (!oBB.isInside(xyz))
                        continue; // don't search outside active voxel area
                    
                    *pbOutsideActiveBounds = false; // we are still inside the bounds
                    
                    int dx = xyz.x() - xyzCenter.X;
                    int dy = xyz.y() - xyzCenter.Y;
                    int dz = xyz.z() - xyzCenter.Z;
                    if (dx * dx + dy * dy + dz * dz <= r2)
                    {
                        bool bInside = oAccess.getValue(xyz) <= 0.0f;
                        if (bInside != bReferenceInside)
                        {
                            pxyzPoint->X = xyz.x();
                            pxyzPoint->Y = xyz.y();
                            pxyzPoint->Z = xyz.z();
                            return true;
                        }
                    }
                }
            }
        }
        
        return false;
    }
};

} // PicoGK namespace

#endif // PICOGKVDBVOXELS_H_
