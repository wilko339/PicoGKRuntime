//
// SPDX-License-Identifier: Apache-2.0
//
// PicoGK ("peacock") is a compact software kernel for computational geometry,
// specifically for use in Computational Engineering Models (CEM).
//
// For more information, please visit https://picogk.org
//
// PicoGK is developed and maintained by LEAP 71 - © 2023-2024 by LEAP 71
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

#ifndef PICOGKGLPOLYLINE_H_
#define PICOGKGLPOLYLINE_H_

#include "PicoGKTypes.h"
#include <vector>
#include <memory>
#include <assert.h>

namespace PicoGK
{
class PolyLine
{
public:
    PKSHAREDPTR(PolyLine);
    
    PolyLine(const ColorFloat& clr)
    {
        m_clrLines = clr;
    }
    
    ~PolyLine()
    {
    }
    
    int32_t nAddVertex(const Vector3& vec)
    {
        m_oVertices.push_back(vec);
        return int32_t(m_oVertices.size() - 1);
    }

    void GetVertex( int32_t nIndex,
                    Vector3* pvec) const
    {
        assert(nIndex < m_oVertices.size());
        *pvec = m_oVertices[nIndex];
    }

    int32_t nVertexCount() const
    {
        return (int32_t) m_oVertices.size();
    }
    
    void* pVertexData() const
    {
        return (void*) m_oVertices.data();
    }

    ColorFloat clrLines() const
    {
        return m_clrLines;
    }
    
protected:
    std::vector<Vector3>    m_oVertices;
    ColorFloat              m_clrLines;
};
}

#endif
