/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

import QtQuick 2.15
import QtQuick3D 1.15
import QtQuick3D.Materials 1.15

CustomMaterial {
    property bool uEnvironmentMappingEnabled: true
    property bool uShadowMappingEnabled: false
    property vector3d tiling: Qt.vector3d(3, 3, 3)
    property real brushing_strength: 0.5
    property real reflection_stretch: 0.5
    property vector3d metal_color: Qt.vector3d(0.95, 0.95, 0.95)
    property real bump_amount: 0.4

//                                         <Property formalName="Roughness Map U" name="roughness_texture_u" description="Horizontal roughness texture" type="Texture" filter="linear" minfilter="linearMipmapLinear" clamp="repeat" usage="roughness" default="./maps/materials/brushed_full_contrast.png" category="Material"/>
//                                         <Property formalName="Roughness Map V" name="roughness_texture_v" description="Vertical roughness texture" type="Texture" filter="linear" minfilter="linearMipmapLinear" clamp="repeat" usage="roughness" default="./maps/materials/brushed_full_contrast.png" category="Material"/>


    shaderInfo: ShaderInfo {
        version: "330"
        type: "GLSL"
        shaderKey: ShaderInfo.Glossy
    }


    property TextureInput uEnvironmentTexture: TextureInput {
            id: uEnvironmentTexture
            enabled: uEnvironmentMappingEnabled
            texture: Texture {
                id: envImage
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/spherical_checker.png"
            }
    }
    property TextureInput uBakedShadowTexture: TextureInput {
            enabled: uShadowMappingEnabled
            texture: Texture {
                id: shadowImage
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/shadow.png"
            }
    }

    property TextureInput brush_texture: TextureInput {
            enabled: true
            texture: Texture {
                id: brushTexture
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/brushed_full_contrast.png"
            }
    }

    property TextureInput roughness_texture_u: TextureInput {
            enabled: true
            texture: Texture {
                id: roughnessUTexture
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/brushed_full_contrast.png"
            }
    }

    property TextureInput roughness_texture_v: TextureInput {
            enabled: true
            texture: Texture {
                id: roughnessVTexture
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/brushed_full_contrast.png"
            }
    }


    property TextureInput bump_texture: TextureInput {
            enabled: true
            texture: Texture {
                id: bumpTexture
                tilingModeHorizontal: Texture.Repeat
                tilingModeVertical: Texture.Repeat
                source: "maps/brushed_a.png"
            }
    }

    Shader {
        id: aluminumBrushedFragShader
        stage: Shader.Fragment
        shader: "shaders/aluminumBrushed.frag"
    }

    passes: [
        Pass {
            shaders: aluminumBrushedFragShader
        }
    ]
}
