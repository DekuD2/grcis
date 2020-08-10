using FilipRechtorik;
using System.Net;

// Youtube video
//https://youtu.be/fVjdsLvRZic

//////////////////////////////////////////////////
// Rendering params.

Debug.Assert(scene != null);
Debug.Assert(context != null);

// Override image resolution and supersampling.
context[PropertyName.CTX_WIDTH]         = 640; // 1800, 1920
context[PropertyName.CTX_HEIGHT]        = 480; // 1200, 1080
//context[PropertyName.CTX_SUPERSAMPLING] =   4; //  400,   64

// Tooltip (if script uses values from 'param').
context[PropertyName.CTX_TOOLTIP] = "n=<double> (index of refraction)\rschlick=<double> (Schlick coeff)\rmat={mirror|glass|diffuse}}";

// Params dictionary.
Dictionary<string, string> p = Util.ParseKeyValueList(param);

// schlick = <schlick-coeff>
double sch = 1.0;
Util.TryParse(p, "schlick", ref sch);

// n = <index-of-refraction>
double n = 1.6;
Util.TryParse(p, "n", ref n);

// mat = {mirror|glass|diffuse}
PhongMaterial pm = new PhongMaterial(new double[] { 0.0, 0.2, 0.1 }, 0.05, 0.05, 0.1, 128, sch);
pm.n  = n;

if (p.TryGetValue("mat", out string mat))
  switch (mat)
  {
    case "diffuse":
      pm = new PhongMaterial(new double[] {0.1, 0.1, 0.6}, 0.1, 0.8, 0.2, 16, sch);
      break;

    case "mirror":
      pm = new PhongMaterial(new double[] {1.0, 1.0, 0.8}, 0.01, 0.0, 1.1, 8192, sch);
      break;

    default:
    case "glass":
      break;
  }
    
//////////////////////////////////////////////////
// CSG scene.

AnimatedCSGInnerNode root = new AnimatedCSGInnerNode(SetOperation.Union);
//CSGInnerNode root = new CSGInnerNode(SetOperation.Union);
root.SetAttribute(PropertyName.REFLECTANCE_MODEL, new PhongModel());
root.SetAttribute(PropertyName.MATERIAL, new PhongMaterial(new double[] {1.0, 0.7, 0.1}, 0.3, 0.4, 0.3, 128, sch));
scene.Intersectable = root;

// Background color.
scene.BackgroundColor = new double[] {0.0, 0.01, 0.03};
scene.Background = new DefaultBackground(scene.BackgroundColor);

// Camera.
scene.Camera = new StaticCamera(new Vector3d(0.7, 0.5, -6.0),
                                new Vector3d(0.0, -0.18, 1.0),
                                50.0);

// Light sources.
scene.Sources = new System.Collections.Generic.LinkedList<ILightSource>();
scene.Sources.Add(new AmbientLightSource(0.8));
scene.Sources.Add(new PointLightSource(new Vector3d(-5.0, 4.0, -3.0), 1.2));

//AnimatedRayScene ascene = scene as AnimatedRayScene;
//if (ascene != null)
//  ascene.End = 20.0;

// --- NODE DEFINITIONS ----------------------------------------------------

// Transparent/mirror/diffuse sphere.

int x = 50;
int y = 100;
int z = 50;

// Create the fuel for our fire
Volume fuel = new Volume(x, y, z);
fuel.SetSeed(10); // This is necessary to do because of a multi-threading bug (idk why but it doesn't always clone, sometimes it creates completely new instances and clones those)
fuel.FuelRandomize(5, y / 4);
// Create fire from our fuel
Fire f = new Fire(x, y, z, fuel);
// And join them in a fire simulation. This is a necessary step for fires that have fuel.
// For fires that don't have fuel and instead have Fire.FuelType == { Floor | Point },
// you can just insert the fire directly into the root node.
FireSimulation sim = new FireSimulation(f, fuel);

pm = new PhongMaterial(new double[] {0.1, 0.1, 0.6}, 0.5, 0.5, 0.3, 16, sch){Kt = 0.1};
sim.SetAttribute(PropertyName.MATERIAL, pm);
root.InsertChild(sim, Matrix4d.RotateY(0) * Matrix4d.Scale(2, 4, 2) * Matrix4d.CreateTranslation(0, 0.6, -0.5));

// Infinite plane with checker.
Plane pl = new Plane();
pl.SetAttribute(PropertyName.COLOR, new double[] {0.3, 0.0, 0.0});
pl.SetAttribute(PropertyName.TEXTURE, new CheckerTexture(0.6, 0.6, new double[] {1.0, 1.0, 1.0}));
root.InsertChild(pl, Matrix4d.RotateX(-MathHelper.PiOver2) * Matrix4d.CreateTranslation(0.0, -1.0, 0.0));
