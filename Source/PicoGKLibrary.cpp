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

#include "PicoGKTypes.h"
#include "PicoGK.h"

#include "PicoGKLibraryMgr.h"
//#include "PicoGKGLViewer.h"

using namespace PicoGK;

void SafeCopyInfoString(const std::string s, char psz[PKINFOSTRINGLEN])
{
#ifdef _WINDOWS
    strncpy_s(psz, PKINFOSTRINGLEN-1, s.c_str(), s.length());
#else
    strncpy(psz, s.c_str(), PKINFOSTRINGLEN-1);
#endif
    psz[PKINFOSTRINGLEN-1] = 0;
}

PICOGK_API void Library_Init(float fVoxelSizeMM, bool bTriangulateMeshes, float fMeshAdaptivity)
{
    Library::oLib().InitLibrary(fVoxelSizeMM, bTriangulateMeshes, fMeshAdaptivity);
}

PICOGK_API void Library_GetName(char psz[PKINFOSTRINGLEN])
{
   SafeCopyInfoString(Library::oLib().strName(), psz);
}

PICOGK_API void Library_GetVersion(char psz[PKINFOSTRINGLEN])
{
    SafeCopyInfoString(Library::oLib().strVersion(), psz);
}

PICOGK_API void Library_GetBuildInfo(char psz[PKINFOSTRINGLEN])
{
    SafeCopyInfoString(Library::oLib().strBuildInfo(), psz);
}

PICOGK_API PKMESH Mesh_hCreate()
{
    return (PKMESH) Library::oLib().proMeshCreate();
}

PICOGK_API PKMESH Mesh_hCreateFromVoxels(PKVOXELS hVoxels)
{
    Voxels::Ptr* proVoxels = (Voxels::Ptr*) hVoxels;
    return (PKMESH) Library::oLib().proMeshCreateFromVoxels(**proVoxels);
}

PICOGK_API bool Mesh_bIsValid(PKMESH hThis)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    return Library::oLib().bMeshIsValid(proThis);
}

PICOGK_API void Mesh_Destroy(PKMESH hThis)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));

    Library::oLib().MeshDestroy(proThis);
}

PICOGK_API int32_t Mesh_nAddVertex( PKMESH hThis,
                                    const Vector3* pvecVertex)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    return (*proThis)->nAddVertex(*pvecVertex);
}
    
PICOGK_API void Mesh_GetVertex( PKMESH      hThis,
                                int32_t     nVertex,
                                Vector3*    pvecVertex)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    (*proThis)->GetVertex(nVertex, pvecVertex);
}

PICOGK_API int32_t Mesh_nVertexCount(PKMESH hThis)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    return (*proThis)->nVertexCount();
}

PICOGK_API int32_t Mesh_nAddTriangle(   PKMESH hThis,
                                        const Triangle* psTri)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    return (*proThis)->nAddTriangle(*psTri);
}

PICOGK_API int32_t Mesh_nAddQuad(       PKMESH hThis,
                                        const Quad* psQuad)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*)hThis;
    assert(Library::oLib().bMeshIsValid(proThis));

    return (*proThis)->nAddQuad(*psQuad);
}

PICOGK_API void Mesh_GetTriangle(   PKMESH hThis,
                                    int32_t nTriangle,
                                    Triangle* psTri)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    return (*proThis)->GetTriangle( nTriangle,
                                    psTri);
}

PICOGK_API void Mesh_GetQuad(       PKMESH hThis,
                                    int32_t nQuad,
                                    Quad* psQuad)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*)hThis;
    assert(Library::oLib().bMeshIsValid(proThis));

    return (*proThis)->GetQuad(nQuad, psQuad);
}

PICOGK_API void Mesh_GetTriangleV(  PKMESH      hThis,
                                    int32_t     nTriangle,
                                    Vector3*    pvecA,
                                    Vector3*    pvecB,
                                    Vector3*    pvecC)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    (*proThis)->GetTriangle(    nTriangle,
                                pvecA,
                                pvecB,
                                pvecC);
}

PICOGK_API void Mesh_GetQuadV(PKMESH hTHis,
    int32_t nQuad,
    PKVector3* pvecA,
    PKVector3* pvecB,
    PKVector3* pvecC,
    PKVector3* pvecD)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*)hTHis;
    assert(Library::oLib().bMeshIsValid(proThis));

    (*proThis)->GetQuad(nQuad,
        pvecA,
        pvecB,
        pvecC,
        pvecD);
}

PICOGK_API void Mesh_GetBoundingBox(    PKMESH hThis,
                                        BBox3* poBox)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    (*proThis)->GetBoundingBox(poBox);
}

PICOGK_API int32_t Mesh_nTriangleCount(PKMESH hThis)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*) hThis;
    assert(Library::oLib().bMeshIsValid(proThis));
    
    return (*proThis)->nTriangleCount();
}

PICOGK_API int32_t Mesh_nQuadCount(PKMESH hThis)
{
    Mesh::Ptr* proThis = (Mesh::Ptr*)hThis;
    assert(Library::oLib().bMeshIsValid(proThis));

    return (*proThis)->nQuadCount();
}

PICOGK_API PKLATTICE Lattice_hCreate()
{
    return (PKLATTICE) Library::oLib().proLatticeCreate();
}

PICOGK_API bool Lattice_bIsValid(PKLATTICE hThis)
{
    Lattice::Ptr* proThis = (Lattice::Ptr*) hThis;
    return Library::oLib().bLatticeIsValid(proThis);
}

PICOGK_API void Lattice_Destroy(PKLATTICE hThis)
{
    Lattice::Ptr* proThis = (Lattice::Ptr*) hThis;
    Library::oLib().LatticeDestroy(proThis);
}

PICOGK_API void Lattice_AddSphere(  PKLATTICE hThis,
                                    const Vector3* vecCenter,
                                    float fRadius)
{
    Lattice::Ptr* proThis = (Lattice::Ptr*) hThis;
    assert(Library::oLib().bLatticeIsValid(proThis));
    
    (*proThis)->AddSphere(  *vecCenter,
                            fRadius);
}

PICOGK_API void Lattice_AddBeam(    PKLATTICE hThis,
                                    const Vector3* pvecA,
                                    const Vector3* pvecB,
                                    float fRadiusA,
                                    float fRadiusB,
                                    bool  bRoundCap)
{
    Lattice::Ptr* proThis = (Lattice::Ptr*) hThis;
    assert(Library::oLib().bLatticeIsValid(proThis));
    
    (*proThis)->AddBeam(    *pvecA,
                            *pvecB,
                            fRadiusA,
                            fRadiusB,
                            bRoundCap);
}

PICOGK_API PKVOXELS Voxels_hCreate()
{
    return (PKVOXELS) Library::oLib().proVoxelsCreate();
}

PICOGK_API PKVOXELS Voxels_hCreateCopy(PKVOXELS hSource)
{
    Voxels::Ptr* proSource = (Voxels::Ptr*) hSource;
    assert(Library::oLib().bVoxelsIsValid(proSource));
    
    return (PKVOXELS) Library::oLib().proVoxelsCreateCopy(**proSource);
}

PICOGK_API bool Voxels_bIsValid(PKVOXELS hThis)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    return Library::oLib().bVoxelsIsValid(proThis);
}

PICOGK_API void Voxels_Destroy(PKVOXELS hThis)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Library::oLib().VoxelsDestroy(proThis);
}

PICOGK_API void Voxels_BoolAdd( PKVOXELS hThis,
                                PKVOXELS hOther)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Voxels::Ptr* proOther = (Voxels::Ptr*) hOther;
    assert(Library::oLib().bVoxelsIsValid(proOther));
    
    (*proThis)->BoolAdd(**proOther);
}

PICOGK_API void Voxels_BoolSubtract( PKVOXELS hThis,
                                     PKVOXELS hOther)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Voxels::Ptr* proOther = (Voxels::Ptr*) hOther;
    assert(Library::oLib().bVoxelsIsValid(proOther));
    
    (*proThis)->BoolSubtract(**proOther);
}

PICOGK_API void Voxels_BoolIntersect(   PKVOXELS hThis,
                                        PKVOXELS hOther)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Voxels::Ptr* proOther = (Voxels::Ptr*) hOther;
    assert(Library::oLib().bVoxelsIsValid(proOther));
    
    (*proThis)->BoolIntersect(**proOther);
}

PICOGK_API void Voxels_Offset(  PKVOXELS hThis,
                                float fDist)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->Offset(fDist, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_DoubleOffset(    PKVOXELS hThis,
                                        float fDist1,
                                        float fDist2)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->DoubleOffset(fDist1, fDist2, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_TripleOffset(    PKVOXELS hThis,
                                        float fDist)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->TripleOffset(fDist, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_RenderMesh(  PKVOXELS hThis,
                                    PKMESH hMesh)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Mesh::Ptr* proMesh = (Mesh::Ptr*) hMesh;
    assert(Library::oLib().bMeshIsValid(proMesh));
    
    (*proThis)->RenderMesh(**proMesh, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_RenderImplicit(  PKVOXELS hThis,
                                        const PKBBox3* poBBox,
                                        PKPFnfSdf pfnSDF)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->RenderImplicit(*poBBox, pfnSDF, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_IntersectImplicit(   PKVOXELS hThis,
                                            PKPFnfSdf pfnSDF)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->IntersectImplicit(pfnSDF, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_RenderLattice(   PKVOXELS hThis,
                                        PKLATTICE hLattice)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Lattice::Ptr* proLattice = (Lattice::Ptr*) hLattice;
    assert(Library::oLib().bLatticeIsValid(proLattice));
    
    (*proThis)->RenderLattice(**proLattice, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_ProjectZSlice( PKVOXELS hThis,
                                      float fZStart,
                                      float fZEnd)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->ProjectZSlice(  fZStart,
                                fZEnd,
                                Library::oLib().fVoxelSizeMM());
}

PICOGK_API bool Voxels_bIsEqual(    PKVOXELS hThis,
                                    PKVOXELS hOther)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    Voxels::Ptr* proOther = (Voxels::Ptr*) hOther;
    assert(Library::oLib().bVoxelsIsValid(proOther));
    
    return (*proThis)->bIsEqual(**proOther);
}

PICOGK_API void Voxels_CalculateProperties( PKVOXELS hThis,
                                            float* pfVolume,
                                            BBox3* poBBox)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->CalculateProperties(pfVolume, poBBox, Library::oLib().fVoxelSizeMM());
}

PICOGK_API void Voxels_GetSurfaceNormal(    PKVOXELS            hThis,
                                            const PKVector3*    pvecSurfacePoint,
                                            PKVector3*          pvecNormal)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    (*proThis)->GetSurfaceNormal(*pvecSurfacePoint, Library::oLib().fVoxelSizeMM(), pvecNormal);
}

PICOGK_API bool Voxels_bClosestPointOnSurface(  PKVOXELS            hThis,
                                                const PKVector3*    pvecSearch,
                                                PKVector3*          pvecSurfacePoint)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    return (*proThis)->bFindClosestPointOnSurface(  *pvecSearch,
                                                    Library::oLib().fVoxelSizeMM(),
                                                    pvecSurfacePoint);
}

PICOGK_API bool Voxels_bRayCastToSurface(   PKVOXELS            hThis,
                                            const PKVector3*    pvecSearch,
                                            const PKVector3*    pvecDirection,
                                            PKVector3*          pvecSurfacePoint)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    return (*proThis)->bRayCastToSurface(   *pvecSearch,
                                            *pvecDirection,
                                            Library::oLib().fVoxelSizeMM(),
                                            pvecSurfacePoint);
}

PICOGK_API void Voxels_GetVoxelDimensions(  PKVOXELS hThis,
                                            int32_t* pnXSize,
                                            int32_t* pnYSize,
                                            int32_t* pnZSize)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    return (*proThis)->GetVoxelDimensions(pnXSize, pnYSize, pnZSize);
}

PICOGK_API void Voxels_GetSlice(    PKVOXELS    hThis,
                                    int32_t     nZSlice,
                                    float*      pfBuffer)
{
    Voxels::Ptr* proThis = (Voxels::Ptr*) hThis;
    assert(Library::oLib().bVoxelsIsValid(proThis));
    
    return (*proThis)->GetSlice(nZSlice, pfBuffer);
}

//PICOGK_API PKVIEWER Viewer_hCreate( const char*             pszWindowTitle,
//                                    const Vector2*          pvecSize,
//                                    PKFInfo                 pfnInfoCallback,
//                                    PKPFUpdateRequested     pfnUpdateCallback,
//                                    PKPFKeyPressed          pfnKeyPressedCallback,
//                                    PKPFMouseMoved          pfnMouseMoveCallback,
//                                    PKPFMouseButton         pfnMouseButtonCallback,
//                                    PKPFScrollWheel         pfnScrollWheelCallback,
//                                    PKPFWindowSize          pfnWindowSize)
//{
//    return (PKVIEWER) ViewerManager::oMgr().poCreate(
//                pszWindowTitle,
//                *pvecSize,
//                pfnInfoCallback,
//                pfnUpdateCallback,
//                pfnKeyPressedCallback,
//                pfnMouseMoveCallback,
//                pfnMouseButtonCallback,
//                pfnScrollWheelCallback,
//                pfnWindowSize);
//}
//
//PICOGK_API bool Viewer_bIsValid(PKVIEWER hThis)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    return ViewerManager::oMgr().bIsValid(poThis);
//}
//
//PICOGK_API void Viewer_Destroy(PKVIEWER hThis)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    PicoGK::ViewerManager::oMgr().Destroy(poThis);
//}
//
//PICOGK_API void Viewer_RequestUpdate(PKVIEWER hThis)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->RequestUpdate();
//}
//
//PICOGK_API bool Viewer_bPoll(PKVIEWER hThis)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    return poThis->bPoll();
//}
//
//PICOGK_API  void Viewer_RequestScreenShot(  PKVIEWER        hThis,
//                                            const char*     pszScreenShotPath)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->RequestScreenShot(pszScreenShotPath);
//}
//
//PICOGK_API void Viewer_RequestClose(PKVIEWER hThis)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->RequestClose();
//}
//
//PICOGK_API bool Viewer_bLoadLightSetup( PKVIEWER        hThis,
//                                        const char*     pDiffTextureDDS,
//                                        int32_t         nDiffTextureSize,
//                                        const char*     pSpecTextureDDS,
//                                        int32_t         nSpecTextureSize)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    return poThis->bLoadLightSetup( pDiffTextureDDS,
//                                    nDiffTextureSize,
//                                    pSpecTextureDDS,
//                                    nSpecTextureSize);
//}
//
//PICOGK_API void Viewer_AddMesh( PKVIEWER    hThis,
//                                int32_t     nGroupID,
//                                PKMESH      hMesh)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    Mesh::Ptr* proMesh = (Mesh::Ptr*) hMesh;
//    assert(Library::oLib().bMeshIsValid(proMesh));
//    
//    poThis->AddMesh(nGroupID, proMesh);
//}
//
//PICOGK_API void Viewer_RemoveMesh(  PKVIEWER hThis,
//                                    PKMESH   hMesh)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    Mesh::Ptr* proMesh = (Mesh::Ptr*) hMesh;
//    assert(Library::oLib().bMeshIsValid(proMesh));
//    
//    poThis->RemoveMesh(proMesh);
//}
//
//PICOGK_API void Viewer_AddPolyLine( PKVIEWER    hThis,
//                                    int32_t     nGroupID,
//                                    PKPOLYLINE  hPolyLine)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    PolyLine::Ptr* proPoly = (PolyLine::Ptr*) hPolyLine;
//    assert(Library::oLib().bPolyLineIsValid(proPoly));
//    
//    poThis->AddPolyLine(    nGroupID,
//                            proPoly);
//}
//
//PICOGK_API void Viewer_RemovePolyLine(  PKVIEWER    hThis,
//                                        PKPOLYLINE  hPolyLine)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    PolyLine::Ptr* proPoly = (PolyLine::Ptr*) hPolyLine;
//    assert(Library::oLib().bPolyLineIsValid(proPoly));
//    
//    poThis->RemovePolyLine(proPoly);
//}
//
//
//PICOGK_API void Viewer_SetGroupVisible( PKVIEWER    hThis,
//                                        int32_t     nGroupID,
//                                        bool        bVisible)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->SetGroupVisible(nGroupID, bVisible);
//}
//
//PICOGK_API void Viewer_SetGroupStatic(  PKVIEWER    hThis,
//                                        int32_t     nGroupID,
//                                        bool        bStatic)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->SetGroupStatic(nGroupID, bStatic);
//}
//
//PICOGK_API void Viewer_SetGroupMaterial(    PKVIEWER            hThis,
//                                            int32_t             nGroupID,
//                                            const ColorFloat*   pclr,
//                                            float               fMetallic,
//                                            float               fRoughness)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->SetGroupMaterial(nGroupID, *pclr, fMetallic, fRoughness);
//}
//
//PICOGK_API void Viewer_SetGroupMatrix(  PKVIEWER            hThis,
//                                        int32_t             nGroupID,
//                                        const Matrix4x4*    pmat)
//{
//    Viewer* poThis = (Viewer*) hThis;
//    assert(ViewerManager::oMgr().bIsValid(poThis));
//    
//    poThis->SetGroupMatrix(nGroupID, *pmat);
//}

PICOGK_API PKVDBFILE VdbFile_hCreate()
{
    return (PKVDBFILE) Library::oLib().proVdbFileCreate();
}

PICOGK_API PKVDBFILE VdbFile_hCreateFromFile(const char* pszFileName)
{
    return (PKVDBFILE) Library::oLib().proVdbFileCreateFromFile(pszFileName);
}

PICOGK_API bool VdbFile_bIsValid(PKVDBFILE hThis)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    return Library::oLib().bVdbFileIsValid(proThis);
}

PICOGK_API void VdbFile_Destroy(PKVDBFILE hThis)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    Library::oLib().VdbFileDestroy(proThis);
}

PICOGK_API bool VdbFile_bSaveToFile(    PKVDBFILE       hThis,
                                        const char*     pszFileName)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    return (*proThis)->bSaveToFile(pszFileName);
}

PICOGK_API PKVOXELS VdbFile_hGetVoxels( PKVDBFILE   hThis,
                                        int32_t     nIndex)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    return (PKVOXELS) Library::oLib().proVdbFileGetVoxels(*proThis, nIndex);
}

PICOGK_API int32_t VdbFile_nAddVoxels(  PKVDBFILE   hThis,
                                        const char* pszFieldName,
                                        PKVOXELS    hVoxels)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    Voxels::Ptr* proVoxels = (Voxels::Ptr*) hVoxels;
    assert(Library::oLib().bVoxelsIsValid(proVoxels));
    
    return Library::oLib().nVdbFileAddVoxels(   *proThis,
                                                pszFieldName,
                                                *proVoxels);
}

PICOGK_API int32_t VdbFile_nFieldCount(PKVDBFILE hThis)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    return (*proThis)->nGridCount();
}

PICOGK_API void VdbFile_GetFieldName(   PKVDBFILE   hThis,
                                        int32_t     nIndex,
                                        char        psz[PKINFOSTRINGLEN])
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    SafeCopyInfoString( (*proThis)->strNameAt(nIndex),
                        psz);
}

PICOGK_API int VdbFile_nFieldType(  PKVDBFILE   hThis,
                                    int32_t     nIndex)
{
    VdbFile::Ptr* proThis = (VdbFile::Ptr*) hThis;
    assert(Library::oLib().bVdbFileIsValid(proThis));
    
    return (*proThis)->nTypeAt(nIndex);
}

