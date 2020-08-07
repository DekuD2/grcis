﻿using MathSupport;
using OpenTK;
using Scene3D;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using Rendering;
using System.Runtime.CompilerServices;
using System.Collections.Specialized;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Threading;
using System.Windows.Forms.VisualStyles;
using System.Security.AccessControl;
using System.CodeDom.Compiler;
using OpenTK.Graphics.ES11;
using System.Security.Cryptography;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Drawing;
using System.Configuration;
using Utilities;
using System.Drawing.Drawing2D;
using OpenTK.Graphics.OpenGL;

namespace FilipRechtorik
{
  public static class PlasmaColor
  {
    // Following links are not used...
    // I just used lerp between red and green in the end...
    // https://en.wikipedia.org/wiki/Color_temperature#Calculation
    // https://en.wikipedia.org/wiki/CIE_1960_color_space
    // https://www.osapublishing.org/DirectPDFAccess/483495D3-02CF-CAB6-FA39786925B7B494_344803/oe-24-13-14066.pdf?da=1&id=344803&seq=0&mobile=no

    static Matrix3d colorTransform = new Matrix3d(
      3.1956, 2.4478, -0.1434,
      -2.5455, 7.0492, 0.9963,
      0.0000, 0.0000, 1.0000);

    static double xe = 0.3320;
    static double ye = 0.1858;

    static Vector3d UVToRGB (double u, double v)
    {
      double x = 3 * u / (2 * u - 8 * v + 4);
      double y = 2 * v / (2 * u - 8 * v + 4);

      return Vector3d.Transform(new Vector3d(x, y, 1), new Matrix4d(colorTransform));
    }

    static double[] flame0 = new double[]{ 89 / 255d, 22 / 255d, 3 / 255d };
    static double f1s = 1;
    static double[] flame1 = new double[]{ f1s * 250 / 255d, f1s * 236 / 255d, f1s * 165 / 255d };

    public static double[] ColorAt (double temperature)
    {
      double Lerp (double from, double to, double k)
        => (1 - k) * from + k * to;

      temperature = temperature * temperature; // * temperature;

      return new double[]
      {
        Lerp(flame0[0], flame1[0], temperature),
        Lerp(flame0[1], flame1[1], temperature),
        Lerp(flame0[2], flame1[2], temperature)
      };
    }
  }

  public class Wind
  {
    Vector3d frequency;
    Vector3d amplitude;
    Vector3d baseShift;
    Vector3d timeShift;

    public Wind ()
    {
      Random rnd = new Random();
      frequency = RandomVector(rnd, new Vector3d(0.3), new Vector3d(0.6));
      amplitude = RandomVector(rnd, new Vector3d(0.03), new Vector3d(0.12));
      baseShift = RandomVector(rnd, new Vector3d(0), new Vector3d(3.14159265358979));
      timeShift = RandomVector(rnd, new Vector3d(-0.05), new Vector3d(0.05));
    }

    Vector3d RandomVector (Random rnd, Vector3d min, Vector3d max)
      => new Vector3d(
        rnd.NextDouble() * (max.X - min.X) + min.X,
        rnd.NextDouble() * (max.Y - min.Y) + min.Y,
        rnd.NextDouble() * (max.Z - min.Z) + min.Z);

    public Vector3d ForceAt (Vector3d position, double time)
    {
      Vector3d wind = new Vector3d(
        Math.Sin(frequency.X * (position.X + baseShift.X + timeShift.X * time)) * amplitude.X,
        Math.Sin(frequency.Y * (position.Y + baseShift.Y + timeShift.Y * time)) * amplitude.Y,
        Math.Sin(frequency.Z * (position.Z + baseShift.Z + timeShift.Z * time)) * amplitude.Z);
      return wind;
    }
  }

  [Serializable]
  public class Fire : Texture3D, ITimeDependent
  {
    public enum Fuel { Floor, Point, Volume };

    public double PlasmaIntensity = 4; //17;
    public double SmokeIntensity = 10;
    public double SamplingFrequency = 10;
    const double smokeExponentialBase = 0.5;

    public Volume FuelVolume;
    public Fuel FuelType = Fuel.Point;

    Wind wind = new Wind ();

    public double AnimationStep = 1 / 25d;
    double simulationTime = 0d;
    public double Start { get; set; }
    public double End { get; set; }
    private double time;
    public double Time
    {
      get => time;
      set
      {
        time = value;
        AnimateTo(time);
      }
    }

    public Fire (int sx, int sy, int sz, Volume fuel = null) : base(sx, sy, sz)
    {
      if (fuel != null)
      {
        FuelType = Fuel.Volume;
        FuelVolume = fuel;
      }
    }

    public void ApplyHeatTransfer (double time)
    {
      void ApplyHeat (float value, float x, float y, float z)
      {
        // Lower coordinates in the texture
        int lx = (int)x;
        int ly = (int)y;
        int lz = (int)z;

        // Uper coordinates
        int ux = lx + 1;
        int uy = ly + 1;
        int uz = lz + 1;

        // Proportions from current to next value
        float px = x - lx;
        float py = y - ly;
        float pz = z - lz;

        // Inverse proportions
        float ipx = 1 - px;
        float ipy = 1 - py;
        float ipz = 1 - pz;

        volume[lx, ly, lz] += ipx * ipy * ipz * value;
        volume[lx, ly, uz] += ipx * ipy * pz * value;
        volume[lx, uy, lz] += ipx * py * ipz * value;
        volume[lx, uy, uz] += ipx * py * pz * value;
        volume[ux, ly, lz] += px * ipy * ipz * value;
        volume[ux, ly, uz] += px * ipy * pz * value;
        volume[ux, uy, lz] += px * py * ipz * value;
        volume[ux, uy, uz] += px * py * pz * value;
      }

      // Create a copy and clear the volume texture
      float[,,] copy = new float[Sx, Sy, Sz];
      for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++)
          for (int k = 0; k < Sz; k++)
          {
            copy[i, j, k] = volume[i, j, k];
            volume[i, j, k] = 0;
          }

      for (int i = 1; i < Sx - 1; i++)
        for (int j = 1; j < Sy - 1; j++)
          for (int k = 1; k < Sz - 1; k++)
          {
            Vector3d f = wind.ForceAt(new Vector3d(i / (Sx - 1), j / (Sy- 1), k / (Sz- 1)), simulationTime);
            f.Y += 6;
            f *= AnimationStep;

            f = Vector3d.Clamp(f, new Vector3d(-0.99), new Vector3d(0.99));
            ApplyHeat(copy[i, j, k], i + (float)f.X, j + (float)f.Y, k + (float)f.Z);
          }
    }

    public void AnimateTo (double time)
    {
      // Idea:
      // Simulate plasma values using the usual grid simulation
      // But simulate smoke differently. Move it chaotically.
      // Possibly even create another level for smoke that will be simulated differently.
      while (simulationTime < time)
      {
        simulationTime += AnimationStep;

        switch (FuelType)
        {
          case Fuel.Floor:
            BurnFloor(1f);
            break;
          case Fuel.Volume:
            BurnFuel(18, 0.5);
            break;
          default:
            BurnPoint(75f);
            break;
        }
        ApplyHeatTransfer(simulationTime);

        float[,,] copy = new float[Sx, Sy, Sz];
        for (int i = 0; i < Sx; i++)
          for (int j = 0; j < Sy; j++)
            for (int k = 0; k < Sz; k++)
              copy[i, j, k] = volume[i, j, k];

        float neighbourScale = 1; //0.5f; // 6 of these
        float corner1Scale = 1; //0.25f; // 12 of these
        float corner2Scale = 1; //0.125f; // 8 of these
        float centerScale = 27f - (6 * neighbourScale + 12 * corner1Scale + 8 * corner2Scale);

        for (int i = 1; i < Sx - 1; i++)
          for (int j = 1; j < Sy - 1; j++)
            for (int k = 1; k < Sz - 1; k++)
            {
              volume[i, j, k] = 0;
              for (int x = -1; x <= 1; x++)
                for (int y = -1; y <= 1; y++)
                  for (int z = -1; z <= 1; z++)
                  {
                    // I want to slow down the dissipation of plasma
                    int distance = Math.Abs(x) + Math.Abs(y) + Math.Abs(z);
                    float scale;
                    switch (distance)
                    {
                      case 0:
                        scale = centerScale;
                        break;
                      case 1:
                        scale = neighbourScale;
                        break;
                      case 2:
                        scale = corner1Scale;
                        break;
                      default:
                        scale = corner2Scale;
                        break;
                    }
                    volume[i, j, k] += copy[i + x, j + y, k + z] * scale;
                  }

              volume[i, j, k] /= 27f;
            }
      }
    }

    void BurnFloor (float value)
    {
      for (int i = 1; i < Sx - 1; i++)
        for (int k = 1; k < Sy - 1; k++)
          volume[i, 1, k] = value;
    }

    void BurnPoint (float value)
    {
      volume[Sx / 2, 1, Sz / 2] = value;
    }

    void BurnFuel (double plasmaPerRadius, double burnPerSecond, double burnThreshold = 0.5)
    {
      // First create plasma from wood
      for (int i = 1; i < (Sx - 1); i++)
        for (int j = 1; j < (Sy - 1); j++)
          for (int k = 1; k < (Sz - 1); k++)
          {
            if (FuelVolume.RadiusAt(i, j, k) > 0 && ( // fuel at this point and air in some direction
                FuelVolume.RadiusAt(i + 1, j, k) == 0 || FuelVolume.RadiusAt(i - 1, j, k) == 0 ||
                FuelVolume.RadiusAt(i, j + 1, k) == 0 || FuelVolume.RadiusAt(i, j - 1, k) == 0 ||
                FuelVolume.RadiusAt(i, j, k + 1) == 0 || FuelVolume.RadiusAt(i, j, k - 1) == 0))
              volume[i, j, k] += (float)(plasmaPerRadius * AnimationStep);
          }

      // Then burn wood by plasma
      for (int i = 1; i < (Sx - 1); i++)
        for (int j = 1; j < (Sy - 1); j++)
          for (int k = 1; k < (Sz - 1); k++)
          {
            var plasma = volume[i,j,k];
            if (plasma > burnThreshold)
              FuelVolume.BurnAt(i, j, k, plasma * burnPerSecond * AnimationStep);
          }
    }

    double FindStart (Vector3d p0, Vector3d p1)
    {
      Vector3d gridPos = new Vector3d((p0.X + 0.5) * (Sx - 1), (p0.Y + 0.5) * (Sy - 1), (p0.Z + 0.5) * (Sz - 1));
      Vector3d gridDir = new Vector3d(p1.X * (Sx - 1), p1.Y * (Sy - 1), p1.Z * (Sz - 1));

      double dx = 0; // distance of x to get into the cube "X" region (between X = -0.5 and X = 0.5 parallel planes)
      double dy = 0;
      double dz = 0;

      // Add some epsilon to get inside the grid and not stay on the edge
      // Before I did this, I would end up on the edge of the cube, sometimes on the wrong side due to numerical imprecision, and that would create a noise effect.
      double e = 0.1;
      if (gridPos.X > (Sx - 1) && gridDir.X < 0)
        dx = (Sx - 1) - gridPos.X - e;
      else if (gridPos.X < 0 && gridDir.X > 0)
        dx = -gridPos.X + e;

      if (gridPos.Y > (Sy - 1) && gridDir.Y < 0)
        dy = (Sy - 1) - gridPos.Y - e;
      else if (gridPos.Y < 0 && gridDir.Y > 0)
        dy = -gridPos.Y + e;

      if (gridPos.Z > (Sz - 1) && gridDir.Z < 0)
        dz = (Sz - 1) - gridPos.Z - e;
      else if (gridPos.Z < 0 && gridDir.Z > 0)
        dz = -gridPos.Z + e;

      double t = dx / gridDir.X; // > 0 since both have same sign
      t = Math.Max(t, dy / gridDir.Y);
      t = Math.Max(t, dz / gridDir.Z);

      return t;
    }

    public override LinkedList<Intersection> Intersect (Vector3d p0, Vector3d p1)
    {
      LinkedList<Intersection> result = new LinkedList<Intersection>();

      double fuelT = -1;
      if (FuelType == Fuel.Volume)
        fuelT = FuelVolume?.Intersect(p0, p1).FirstOrDefault()?.T ?? -1;

      double t = FindStart(p0, p1);

      Vector3d pStart = p0 + t * p1;

      if (t > 0 &&
        pStart.X >= -0.5 && pStart.X <= 0.5 &&
        pStart.Y >= -0.5 && pStart.Y <= 0.5 &&
        pStart.Z >= -0.5 && pStart.Z <= 0.5)
        result.AddLast(new Intersection(this)
        {
          T = t,
          Enter = true,
          Front = true,
          NormalLocal = -p1,
          CoordLocal = pStart,
          TangentU = new Vector3d(fuelT)
        })
      ;

      return result;
    }

    /// <summary>
    /// Complete all relevant items in the given Intersection object.
    /// </summary>
    /// <param name="inter">Intersection instance to complete.</param>
    public override void CompleteIntersection (Intersection inter)
    {
      Vector3d p1 = -inter.NormalLocal;
      Vector3d p0 = inter.CoordLocal;

      double fuelT = inter.TangentU.X;

      // How many lines?
      double step = 1d / (SamplingFrequency * maxDim * p1.Length); // unit one / max dimension
      double stepSize = step * p1.Length;

      bool IsInside (Vector3d point)
      {
        return
          point.X >= -0.5 && point.Y >= -0.5 && point.Z >= -0.5 &&
          point.X <= 0.5 && point.Y <= 0.5 && point.Z <= 0.5;
      }

      // Idea for color calculating
      //  Have a smoke factor starting at 0
      //  When you sample value < plasmaThreshold, add it to smoke factor
      //  When you sample value >= plasmaThreshold, compute it's color from gradient, apply smoke factor, add
      //  Return smoke factor as .TextureCoord.X, return color as .SurfaceColor


      Vector3d p = p0;
      double t = 0;
      double smoke = 0;
      double plasma = 0;
      double plasmaThreshold = 0.5;
      double smokeThreshold = 0.1;

      double plasmaMultiplier = (1d / (1 - plasmaThreshold));
      double smokeMultiplier = (1d / (plasmaThreshold - smokeThreshold));
      double[] plasmaColor = new double[3];

      while (IsInside(p) && (fuelT == -1 || t < fuelT))
      {
        double val = Sample(p.X + 0.5, p.Y + 0.5, p.Z + 0.5);

        if (val > plasmaThreshold)
        {
          plasma += val;
          var col = PlasmaColor.ColorAt((val - plasmaThreshold) * plasmaMultiplier);
          double smokeAbsorb = Math.Pow(smokeExponentialBase, smoke);
          plasmaColor[0] += col[0] * smokeAbsorb * PlasmaIntensity * stepSize;
          plasmaColor[1] += col[1] * smokeAbsorb * PlasmaIntensity * stepSize;
          plasmaColor[2] += col[2] * smokeAbsorb * PlasmaIntensity * stepSize;
        }
        else if (val > smokeThreshold)
        {
          smoke += (val - smokeThreshold) * smokeMultiplier * SmokeIntensity * stepSize;
        }

        t += step;
        p = p0 + p1 * t;
      }

      inter.SurfaceColor = plasmaColor;
      inter.TextureCoord.X = plasma;
      inter.TextureCoord.Y = Math.Pow(smokeExponentialBase, smoke);
    }

    public object Clone ()
    {
      //Debug.WriteLine("Clone {");
      //Debug.Indent();
      Fire copy = new Fire(Sx, Sy, Sz, (Volume)FuelVolume?.Clone());

      ShareCloneAttributes(copy);

      copy.PlasmaIntensity = PlasmaIntensity;
      copy.SmokeIntensity = SmokeIntensity;
      copy.SamplingFrequency = SamplingFrequency;
      copy.AnimationStep = AnimationStep;
      copy.simulationTime = simulationTime;
      copy.wind = wind;
      copy.FuelType = FuelType;
      //copy.Fuel = (Volume)Fuel.Clone();
      //copy.Fuel = Fuel;

      for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++)
          for (int k = 0; k < Sz; k++)
            copy.volume[i, j, k] = volume[i, j, k];

      copy.time = time;
      copy.seed = seed;
      copy.random = new Random(seed);

      //Debug.Unindent();
      //Debug.WriteLine("}");
      return copy;
    }

    public static long RecursionFunction (Intersection i, Vector3d dir, double importance, out RayRecursion rr)
    {
      // What if I literally insert the wood inside the fire and then let a ray hit it?
      double plasma = i.TextureCoord.X / 100d;
      double ip = MathHelper.Clamp(1.4 - plasma, 0, 1);
      double smoke = MathHelper.Clamp(i.TextureCoord.Y, 0, 1);
      double fuelT = i.TangentU.X;

      //plasma = 1d;
      //fuelT = -1;

      double length = Math.Sqrt(i.SurfaceColor[0] * i.SurfaceColor[0] + i.SurfaceColor[1] * i.SurfaceColor[1] + i.SurfaceColor[2] * i.SurfaceColor[2]);
      //double l = length;
      //double l2 = 800 - length;
      //if (l > 0)
      //  Debug.WriteLine(l);
      //length = Math.Max(1, length);
      //length = plasma;
      //length = 1;
      double[] col = new double[] { i.SurfaceColor[0] / length, i.SurfaceColor[1] / length, i.SurfaceColor[2] / length };
      col = new double[] { col[0] * plasma, col[1] * plasma, col[2] * plasma };
      if (plasma == 0)
        col = new double[] { 0, 0, 0 };

      if (fuelT == -1)
      {
        rr = new RayRecursion(
          //Util.ColorClone(col, 0.8),
          Util.ColorClone(col),
          //Util.ColorClone(i.SurfaceColor),
          new RayRecursion.RayContribution(i, dir, importance)
          {
            coefficient = new double[] { smoke, smoke, smoke }
          });
        //if (smoke < 0.7)
        //  Debug.WriteLine(smoke);
      }
      else
      {
        double[] wood = new double[] { 59 / 255d, 37 / 255d, 14 / 255d };
        //Debug.WriteLine(l2);
        //rr = new RayRecursion(Util.ColorClone(col, 0.8));
        //rr = new RayRecursion(new double[] { i.SurfaceColor[0] + wood[0], i.SurfaceColor[1] + wood[1], i.SurfaceColor[2] + wood[2] });
        rr = new RayRecursion(new double[] { col[0] + wood[0] * ip, col[1] + wood[1] * ip, col[2] + wood[2] * ip });
      }


      return 122L;
    }
  }
}
