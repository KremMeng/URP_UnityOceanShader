using UnityEngine;
using UnityEngine.Rendering.Universal;

// [ExecuteAlways]
public class PlanerReflection : MonoBehaviour
{
        [SerializeField] public int textureSize = 1024;

        [SerializeField] private Transform ReflectionPlane;
        [SerializeField] public Vector3 planeNormal => ReflectionPlane.up;  
        [SerializeField] public Vector3 planePosition => ReflectionPlane.position; //?
    

        private RenderTexture _reflectionRT;
        
        private Camera _reflectionCamera; 
        private static readonly int ReflectionTexture = Shader.PropertyToID("_ReflectionTexture");


        //creating the reflection camera and RT
        private void Start()
        {
            CreateReflectCamera();
        }
        
        //call before cams start to render,update would late 1 fame
        private void OnPreRender()
        {
            
            UpdateReflectCamParams(Camera.current, _reflectionCamera);
            RenderReflectCam();
        }
        
        private void CreateReflectCamera()
        {
            if (_reflectionCamera == null && _reflectionRT == null)
            {
                
                GameObject gameobject = new GameObject("Reflection Camera");
                _reflectionCamera = gameobject.AddComponent<Camera>();
                _reflectionCamera.CopyFrom(Camera.current);
                
                _reflectionRT = RenderTexture.GetTemporary(textureSize, textureSize, 16, RenderTextureFormat.R8);
                _reflectionCamera.targetTexture = _reflectionRT; //render to RT

                _reflectionCamera.enabled = false;
                Shader.SetGlobalTexture(ReflectionTexture,_reflectionRT);
            }
        }

        private void UpdateReflectCamParams(Camera mainCamera,Camera reflectionCamera)
        {
            if (mainCamera == null && reflectionCamera == null) return;
            reflectionCamera.CopyFrom(mainCamera);
            reflectionCamera.useOcclusionCulling = false;
            if (reflectionCamera.gameObject.TryGetComponent(out UniversalAdditionalCameraData cameraData))
            {
                cameraData.renderShadows = false;
            }
        }
        
        private void RenderReflectCam()
        {
            //w2c * reflectionMatrix
            var reflectionMatrix = ReflectionMatrix(planePosition, planeNormal);
            _reflectionCamera.worldToCameraMatrix = Camera.current.worldToCameraMatrix * reflectionMatrix;
            
            //Projection and Cull underwater View
            CalculateClipPlane();
            
            //Render the camera
            //Matrix reverse Vertexes only,not include Normals
            GL.invertCulling = true; 
            _reflectionCamera.Render();
            GL.invertCulling = false;
        }

        private void CalculateClipPlane()
        {
            var d = -Vector3.Dot(planeNormal, planePosition);
            var plane =new Vector4(planeNormal.x, planeNormal.y, planeNormal.z, d);
            _reflectionCamera.projectionMatrix = Camera.current.CalculateObliqueMatrix(plane);
        }
        private Matrix4x4 ReflectionMatrix(Vector3 orign,Vector3 normal)
        {
            Matrix4x4 reflectionMatrix = Matrix4x4.identity;
            float d = -Vector3.Dot(normal, orign);

            reflectionMatrix.m00 = 1 - 2 * normal.x * normal.x;
            reflectionMatrix.m01 = -2 * normal.x * normal.y;
            reflectionMatrix.m02 = -2 * normal.x * normal.z;
            reflectionMatrix.m03 = -2 * d * normal.x;
            reflectionMatrix.m10 = -2 * normal.x * normal.y;
            reflectionMatrix.m11 = 1 - 2 * normal.y * normal.y;
            reflectionMatrix.m12 = -2 * normal.y * normal.z;
            reflectionMatrix.m13 = -2 * d * normal.y;
            reflectionMatrix.m20 = -2 * normal.x * normal.z;
            reflectionMatrix.m21 = -2 * normal.y * normal.z;
            reflectionMatrix.m22 = 1 - 2 * normal.z * normal.z;
            reflectionMatrix.m23 = -2 * d * normal.z;
            reflectionMatrix.m30 = 0;
            reflectionMatrix.m31 = 0;
            reflectionMatrix.m32 = 0;
            reflectionMatrix.m33 = 1;

            return reflectionMatrix;
        }
    }
