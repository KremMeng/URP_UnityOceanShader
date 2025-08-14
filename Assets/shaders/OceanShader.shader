Shader "Universal Render Pipeline/OceanShader"
{
    Properties
    {
        _MainTex ("mainTex", 2D) = "white" {}
        _ReflectionColor("ReflectionColor",Color) = (1,1,1,1)
        _Diffuse("Diffuse",Color) = (1,1,1,1)
        _Specular("Specular",Color) = (0.77,0.66,0.52,1)
        _Glossness("Glossness",Range(400,800)) = 800
        _NormalBias("NormalBias",Range(0,1)) = 0.7
        _RefactionDistortion("RefractionDistortion",Range(0,1)) = 0.5
        _BackLightTense("BackLightTense",Range(0,500)) = 30
        _EdgeArea("EdgeArea",Range(5,20))=8
        _TransparentFactor("Transparency",Range(0,1)) = 0.28
        _DeepColor("DeepColor",Color) = (0.1, 0.5, 0.7, 1)
        _ShallowColor("ShallowColor",Color)=(0.275, 0.855, 1, 1)
    }
    SubShader
    {
        Tags {
                "RenderPipeline" = "UniversalPipeline"
                "RenderType" = "Transparent"
               // "IgnoreProjector" = "True"
                "Queue" = "Transparent"
        }
        LOD 100
        Blend SrcAlpha OneMinusSrcAlpha
        ZWrite Off

        HLSLINCLUDE
         #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
         #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
         #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Input.hlsl"
         #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/DeclareDepthTexture.hlsl"
         #include "Packages/com.unity.render-pipelines.core/ShaderLibrary/SpaceTransforms.hlsl"
         #include "Packages/com.unity.render-pipelines.core/ShaderLibrary/Common.hlsl"

         //sampler2D _MainTex;
         TEXTURE2D(_MainTex);
         SAMPLER(sampler_MainTex);
         
         //TEXTURE2D(_CameraDepthTexture);
         //SAMPLER(sampler_cameraDepthTexture);
         TEXTURE2D(_CameraOpaqueTexture);
         SAMPLER(sampler_CameraOpaqueTexture);
         
         TEXTURE2D(_depthLutTex);
         SAMPLER(sampler_depthLutTex);

         TEXTURE2D(_ReflectionTexture);
         SAMPLER(sampler_ReflectionTexture);
         
        //暴露在编辑器上的变量需要缓存到常量buffer
         CBUFFER_START(UnityPerMaterial)
            half4 _Diffuse;
            half4 _Specular;
            float _Glossness;
            float4 _MainTex_ST;
            float _NormalBias;
            float _RefactionDistortion;
            float _BackLightTense;
            float _EdgeArea;
            float _TransparentFactor;
            float3 _DeepColor;
            float3 _ShallowColor;
            float4 _Reflection;
         CBUFFER_END
        ENDHLSL
        
        Pass
        {
            //Tags {"LightMode" = "UniversalForward"}
            HLSLPROGRAM
            #pragma vertex Vertex
            #pragma fragment Pixel
            //#pragma multi_compile_fog
            #pragma multi_compile_instancing
            
            // make fog work
            #pragma multi_compile_fog

            struct VertexInput
            {
                float4 vertex : POSITION;
                float3 normal : NORMAL;
            };

            struct VertexOutput
            {
                float4 pos : SV_POSITION;
                float3 worldNormal : TEXCOORD0;
                float3 worldPos : TEXCOORD1;
                float4 screenPos : TEXCOORD2;
                float3 color : COLOR;
                float2 uv : TEXCOORD3;
                float3 fogCoord :TEXCOORD4;
            };
            
            
            VertexOutput Vertex(VertexInput v)
            {
                VertexOutput o = (VertexOutput)0;
                //VertexOutput o;
          
                o.color = _Diffuse.rgb;
                o.pos = TransformObjectToHClip(v.vertex.xyz);
                o.screenPos = ComputeScreenPos(o.pos);
                o.worldNormal = TransformObjectToWorldNormal(v.normal);
                o.worldPos = TransformObjectToWorld(v.vertex.xyz);
                o.fogCoord = o.pos.xyz;
                return o;
            }
            

            half4 Pixel(VertexOutput i) : SV_Target
            {
                //fixed3 ambient = UNITY_LIGHTMODEL_AMBIENT.rgb;
                half3 worldNormal = normalize(i.worldNormal);
                half3 worldLightDir = normalize(_MainLightPosition.rgb);
                //diffuse
                half3 halfLambert = dot(worldNormal,worldLightDir) * 0.5 + 0.5;
                half3 diffuse = _MainLightColor.rgb * _Diffuse.rgb * halfLambert;
                
                //Blinn-Phong
                half3 viewDir = normalize(_WorldSpaceCameraPos.rgb - i.worldPos.rgb);
                half3 H = normalize(worldLightDir + viewDir);
                half3 specular = _MainLightColor.rgb * _Specular.rgb * pow(max(0,dot(worldNormal,H)),_Glossness);

                //Cook-Torrance highlight
                
                
                //transfluency SSS
                half3 H0 = normalize(worldLightDir + _NormalBias * worldNormal);
                half3 backDot = dot(viewDir,-H0);
                half3 backDir = saturate(backDot);
                half3 sssColor = pow(backDir,_EdgeArea) * _BackLightTense;

                float3 baseColor =  diffuse + specular + sssColor ;
                
                //depth based LUT
                float2 screenUV = i.screenPos.xy / i.screenPos.w;
                float depth = SampleSceneDepth(screenUV);
                float viewDepth = LinearEyeDepth(depth,_ZBufferParams); //clip2view
                float linearDepth =  Linear01Depth(viewDepth,_ZBufferParams); //[0,1]
                // float2 lutUV = float2(linearDepth,0.5);
                // float3 lutColor = SAMPLE_TEXTURE2D(_depthLutTex,sampler_depthLutTex,lutUV);
                //float3 lutColor = lerp(_ShallowColor,_DeepColor,linearDepth);
                float alpha = 1.0 - saturate(linearDepth * _TransparentFactor);

                //reflection
                float3 reflectionColor = SAMPLE_TEXTURE2D(_ReflectionTexture,sampler_ReflectionTexture,screenUV).r;

                //refraction
                float2 offset = worldNormal.xy * _RefactionDistortion  * linearDepth * 0.1; //normal is vec3
                screenUV.xy += offset;
                float3 refraction = SAMPLE_TEXTURE2D(_CameraOpaqueTexture,sampler_CameraOpaqueTexture,screenUV).rgb;
                
                //fresnel
                float fresnel = pow(1.0 - max(0, dot(worldNormal, viewDir)), 5.0);
                fresnel = 0.98 * fresnel * _Specular.rgb;
                
                // //absorb && scatter
                // float3 absorption = lerp(_ShallowColor,_DeepColor,linearDepth);
                // float3 scattering = lerp(_DeepColor,_ShallowColor,linearDepth);
                // float3 transfuencyColor = lerp(absorption,scattering,fresnel);
                
                //foam sampling
                
                
                float3 refcol = lerp(refraction,baseColor,saturate(depth * _TransparentFactor));

                float3 finalColor = lerp(refcol,reflectionColor,fresnel);
                
                return half4(finalColor,alpha);
            }
           ENDHLSL
        }
    }
}
