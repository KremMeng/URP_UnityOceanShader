using UnityEngine;

public class CameraDepth : MonoBehaviour
{
    // Start is called before the first frame update
    private void Start()
    {
        Camera cam = GetComponent<Camera>();
        if (cam != null)
        {
            cam.depthTextureMode = DepthTextureMode.Depth;
        }
    }
}
