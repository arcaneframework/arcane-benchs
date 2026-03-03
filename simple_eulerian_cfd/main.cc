// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* main.cc                                                     (C) 2000-2026 */
/*                                                                           */
/* Simple 2D Eulerian CFD.                                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// The following code is a C++ port from the following page:
//
// https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html
//
// with the following Licence:

// Copyright 2022 Matthias Müller - Ten Minute Physics,
// www.youtube.com/c/TenMinutePhysics
// www.matthiasMueller.info/tenMinutePhysics

// MIT License

// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/launcher/ArcaneLauncher.h>

#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/PlatformUtils.h>
#include <arcane/utils/FatalErrorException.h>

#include <arcane/core/MathUtils.h>
#include <arcane/core/ISubDomain.h>

#include <iostream>

using namespace Arcane;

const Int32 U_FIELD = 0;
const Int32 V_FIELD = 1;
const Int32 S_FIELD = 2;
Int32 cnt = 0;

Int32 canvas_width = 600;
Int32 canvas_height = 400;

Real simHeight = 1.1;
Real cScale = canvas_height / simHeight;
Real simWidth = canvas_width / cScale;

class Fluid;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Utils
{
 public:

  static Int32 _doInt32Floor(Real v)
  {
    return static_cast<Int32>(std::floor(v));
  }
  static void _allocateAndFillZero(Real** ptr, Int32 size)
  {
    *ptr = new Real[size];
    Real* p = *ptr;
    for (Int32 i = 0; i < size; ++i)
      p[i] = {};
  }
  static void _allocateAndFillZero(NumArray<Real, MDDim1>* ptr, Int32 size)
  {
    ptr->resize(size);
    SmallSpan<Real> p = ptr->to1DSmallSpan();
    for (Int32 i = 0; i < size; ++i)
      p[i] = {};
  }

  static void _doCopy(Real* new_v, Real* v, Int32 size)
  {
    for (Int32 i = 0; i < size; ++i)
      new_v[i] = v[i];
  }
  static void _doCopy(NumArray<Real, MDDim1>& new_v, NumArray<Real, MDDim1>& v, [[maybe_unused]] Int32 size)
  {
    new_v.copy(v);
  }
  static void _doFillZero(Real* v, Int32 size)
  {
    for (Int32 i = 0; i < size; ++i)
      v[i] = {};
  }
  static void _doFillZero(NumArray<Real, MDDim1>& v, [[maybe_unused]] Int32 size)
  {
    v.fill(0.0);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Scene
{
 public:

  Real gravity = -9.81;
  Real dt = 1.0 / 120.0;
  Int32 numIters = 100;
  Int32 frameNr = 0;
  Real overRelaxation = 1.9;
  Real obstacleX = 0.0;
  Real obstacleY = 0.0;
  Real obstacleRadius = 0.15;
  bool paused = false;
  Int32 sceneNr = 0;
  bool showObstacle = false;
  bool showStreamlines = false;
  bool showVelocities = false;
  bool showPressure = false;
  bool showSmoke = true;
  Fluid* fluid = nullptr;

 public:

  void setupScene(Int32 resolution);
  void setObstacle(Real x, Real y, bool reset);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class Fluid
{
  //using RealArrayType = NumArray<Real, MDDim1>;
  using RealArrayType = Real*;

 public:

  Fluid(Real density_, Int32 numX_, Int32 numY_, Real h_, Scene* scene_)
  {
    density = density_;
    numX = numX_ + 2;
    numY = numY_ + 2;
    h = h_;
    scene = scene_;

    numCells = numX * numY;
    Utils::_allocateAndFillZero(&u, numCells);
    Utils::_allocateAndFillZero(&v, numCells);
    Utils::_allocateAndFillZero(&newU, numCells);
    Utils::_allocateAndFillZero(&newV, numCells);
    Utils::_allocateAndFillZero(&p, numCells);
    Utils::_allocateAndFillZero(&s, numCells);
    Utils::_allocateAndFillZero(&m, numCells);
    Utils::_allocateAndFillZero(&newM, numCells);
  }

  void integrate(Real dt, Real gravity);
  void solveIncompressibility(Int32 numIters, Real dt);
  void extrapolate();
  Real sampleField(Real x, Real y, RealArrayType& f, Int32 field);

  Real avgU(Int32 i, Int32 j)
  {
    Int32 n = numY;
    Real lu = (u[i * n + j - 1] + u[i * n + j] + u[(i + 1) * n + j - 1] + u[(i + 1) * n + j]) * 0.25;
    return lu;
  }

  Real avgV(Int32 i, Int32 j)
  {
    Int32 n = numY;
    Real lv = (v[(i - 1) * n + j] + v[i * n + j] + v[(i - 1) * n + j + 1] + v[i * n + j + 1]) * 0.25;
    return lv;
  }

  void advectVel(Real dt);
  void advectSmoke(Real dt);
  void simulate(Real dt, Real gravity, Int32 numIters);

 public:

  Int32 numX = 0;
  Int32 numY = 0;
  Int32 numCells = 0;
  Real density = 0.0;
  Real h = 0.0;
  RealArrayType u;
  RealArrayType v;
  RealArrayType newU;
  RealArrayType newV;
  RealArrayType p;
  RealArrayType s;
  RealArrayType m;
  RealArrayType newM;
  Scene* scene = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
integrate(Real dt, Real gravity)
{
  Int32 n = numY;
  for (Int32 i = 1; i < numX; i++) {
    for (Int32 j = 1; j < numY - 1; j++) {
      if (s[i * n + j] != 0.0 && s[i * n + j - 1] != 0.0)
        v[i * n + j] += gravity * dt;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
solveIncompressibility(Int32 numIters, Real dt)
{
  Int32 n = numY;
  Real cp = density * h / dt;

  for (Int32 iter = 0; iter < numIters; iter++) {

    for (Int32 i = 1; i < numX - 1; i++) {
      for (Int32 j = 1; j < numY - 1; j++) {

        if (s[i * n + j] == 0.0)
          continue;

        Real sx0 = s[(i - 1) * n + j];
        Real sx1 = s[(i + 1) * n + j];
        Real sy0 = s[i * n + j - 1];
        Real sy1 = s[i * n + j + 1];
        Real ls = sx0 + sx1 + sy0 + sy1;

        if (ls == 0.0)
          continue;

        Real div = u[(i + 1) * n + j] - u[i * n + j] + v[i * n + j + 1] - v[i * n + j];
        Real lp = -div / ls;
        lp *= scene->overRelaxation;
        p[i * n + j] += cp * lp;
        u[i * n + j] -= sx0 * lp;
        u[(i + 1) * n + j] += sx1 * lp;
        v[i * n + j] -= sy0 * lp;
        v[i * n + j + 1] += sy1 * lp;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
extrapolate()
{
  Int32 n = numY;
  for (Int32 i = 0; i < numX; i++) {
    u[i * n + 0] = u[i * n + 1];
    u[i * n + numY - 1] = u[i * n + numY - 2];
  }
  for (Int32 j = 0; j < numY; j++) {
    v[0 * n + j] = v[1 * n + j];
    v[(numX - 1) * n + j] = v[(numX - 2) * n + j];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real Fluid::
sampleField(Real x, Real y, RealArrayType& f, Int32 field)
{
  Int32 n = numY;
  Real h1 = 1.0 / h;
  Real h2 = 0.5 * h;

  x = std::max(std::min(x, numX * h), h);
  y = std::max(std::min(y, numY * h), h);

  Real dx = 0.0;
  Real dy = 0.0;

  switch (field) {
  case U_FIELD:
    dy = h2;
    break;
  case V_FIELD:
    dx = h2;
    break;
  case S_FIELD:
    dx = h2;
    dy = h2;
    break;
  }

  Int32 x0 = std::min(Utils::_doInt32Floor((x - dx) * h1), numX - 1);
  Real tx = ((x - dx) - x0 * h) * h1;
  Int32 x1 = std::min(x0 + 1, numX - 1);

  Int32 y0 = std::min(Utils::_doInt32Floor((y - dy) * h1), numY - 1);
  Real ty = ((y - dy) - y0 * h) * h1;
  Int32 y1 = std::min(y0 + 1, numY - 1);

  Real sx = 1.0 - tx;
  Real sy = 1.0 - ty;

  Real val = sx * sy * f[x0 * n + y0] +
  tx * sy * f[x1 * n + y0] +
  tx * ty * f[x1 * n + y1] +
  sx * ty * f[x0 * n + y1];

  return val;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
advectVel(Real dt)
{
  Utils::_doCopy(newU, u, numCells);
  Utils::_doCopy(newV, v, numCells);

  Int32 n = numY;
  Real h2 = 0.5 * h;

  for (Int32 i = 1; i < numX; i++) {
    for (Int32 j = 1; j < numY; j++) {

      cnt++;

      // u component
      if (s[i * n + j] != 0.0 && s[(i - 1) * n + j] != 0.0 && j < numY - 1) {
        Real x = i * h;
        Real y = j * h + h2;
        Real lu = u[i * n + j];
        Real lv = avgV(i, j);
        x = x - dt * lu;
        y = y - dt * lv;
        Real new_u = sampleField(x, y, u, U_FIELD);
        newU[i * n + j] = new_u;
      }
      // v component
      if (s[i * n + j] != 0.0 && s[i * n + j - 1] != 0.0 && i < numX - 1) {
        Real x = i * h + h2;
        Real y = j * h;
        Real lu = avgU(i, j);
        Real lv = v[i * n + j];
        x = x - dt * lu;
        y = y - dt * lv;
        Real new_v = sampleField(x, y, v, V_FIELD);
        newV[i * n + j] = new_v;
      }
    }
  }

  Utils::_doCopy(u, newU, numCells);
  Utils::_doCopy(v, newV, numCells);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
advectSmoke(Real dt)
{
  Utils::_doCopy(newM, m, numCells);

  Int32 n = numY;
  Real h2 = 0.5 * h;

  for (Int32 i = 1; i < numX - 1; i++) {
    for (Int32 j = 1; j < numY - 1; j++) {

      if (s[i * n + j] != 0.0) {
        Real lu = (u[i * n + j] + u[(i + 1) * n + j]) * 0.5;
        Real lv = (v[i * n + j] + v[i * n + j + 1]) * 0.5;
        Real x = i * h + h2 - dt * lu;
        Real y = j * h + h2 - dt * lv;

        newM[i * n + j] = sampleField(x, y, m, S_FIELD);
      }
    }
  }
  Utils::_doCopy(m, newM, numCells);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
simulate(Real dt, Real gravity, Int32 numIters)
{
  integrate(dt, gravity);

  Utils::_doFillZero(p, numCells);
  solveIncompressibility(numIters, dt);

  extrapolate();
  advectVel(dt);
  advectSmoke(dt);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Scene::
setupScene(Int32 resolution)
{
  obstacleRadius = 0.15;
  overRelaxation = 1.9;

  dt = 1.0 / 60.0;
  numIters = 40;

  Int32 res = 100;

  Real domainHeight = 1.0;
  Real domainWidth = domainHeight / simHeight * simWidth;
  Real h = domainHeight / res;

  Int32 numX = Utils::_doInt32Floor(domainWidth / h);
  Int32 numY = Utils::_doInt32Floor(domainHeight / h);
  std::cout << "domain_height=" << domainHeight << " width=" << domainWidth << " computed_h=" << h << "\n";

  Real density = 1000.0;

  std::cout << "NumX=" << numX << " NumY=" << numY << " h=" << h << "\n";

  Fluid* f_ptr = new Fluid(density, numX, numY, h, this);
  fluid = f_ptr;
  Fluid& f = *f_ptr;

  Int32 n = f.numY;

  {
    // vortex shedding
    Real inVel = 2.0;
    for (Int32 i = 0; i < f.numX; i++) {
      for (Int32 j = 0; j < f.numY; j++) {
        Real s = 1.0; // fluid
        if (i == 0 || j == 0 || j == f.numY - 1)
          s = 0.0; // solid
        f.s[i * n + j] = s;

        if (i == 1) {
          f.u[i * n + j] = inVel;
        }
      }
    }

    Real pipeH = 0.1 * f.numY;
    Int32 minJ = Utils::_doInt32Floor(0.5 * f.numY - 0.5 * pipeH);
    Int32 maxJ = Utils::_doInt32Floor(0.5 * f.numY + 0.5 * pipeH);

    for (Int32 j = minJ; j < maxJ; j++)
      f.m[j] = 0.0;

    setObstacle(0.4, 0.5, true);

    gravity = 0.0;
    showPressure = false;
    showSmoke = true;
    showStreamlines = false;
    showVelocities = false;

    {
      // Vortex shedding
      dt = 1.0 / 120.0;
      numIters = 100;
      showPressure = true;
    }
  }
}

void Scene::
setObstacle(Real x, Real y, bool reset)
{
  Real vx = 0.0;
  Real vy = 0.0;

  if (!reset) {
    vx = (x - obstacleX) / dt;
    vy = (y - obstacleY) / dt;
  }

  obstacleX = x;
  obstacleY = y;
  Real r = obstacleRadius;
  Fluid& f = *fluid;
  Int32 n = f.numY;

  for (Int32 i = 1; i < f.numX - 2; i++) {
    for (Int32 j = 1; j < f.numY - 2; j++) {

      f.s[i * n + j] = 1.0;

      Real dx = (i + 0.5) * f.h - x;
      Real dy = (j + 0.5) * f.h - y;

      if (dx * dx + dy * dy < r * r) {
        f.s[i * n + j] = 0.0;
        f.m[i * n + j] = 1.0;
        f.u[i * n + j] = vx;
        f.u[(i + 1) * n + j] = vx;
        f.v[i * n + j] = vy;
        f.v[i * n + j + 1] = vy;
      }
    }
  }

  showObstacle = true;
}

void executeCode(ISubDomain* sd)
{
  // This code does not use a mesh

  // Get the instance to manage listing.
  ITraceMng* tr = sd->traceMng();

  // Number of cells by dimension
  Int32 resolution = 100;
  // Number of time step
  Int32 nb_time_step = 100;

  Scene scene;
  scene.setupScene(resolution);
  Real current_time = 0.0;
  Real begin_x = Platform::getRealTime();
  for (Int32 i = 0; i < 100; ++i) {
    Fluid& f = *scene.fluid;
    if (i < 10 || ((i % 5) == 0))
      tr->info() << "DoSimulate nx=" << f.numX << " ny=" << f.numY << " dt=" << scene.dt
                 << " gravity=" << scene.gravity << " frame=" << scene.frameNr;
    f.simulate(scene.dt, scene.gravity, scene.numIters);
    ++scene.frameNr;
    current_time += scene.dt;
  }
  Real end_x = Platform::getRealTime();
  tr->info() << "ComputeTime=" << (end_x - begin_x);

  // Check velocity
  // Only valid if resolution == 100 and nb_time_step == 100
  if (nb_time_step == 100 && resolution == 100) {
    Real velocity_sum = 0.0;
    Fluid& f = *scene.fluid;
    Int32 n = f.numY;
    for (Int32 i = 1; i < f.numX - 2; i++) {
      for (Int32 j = 1; j < f.numY - 2; j++) {
        Real u = f.u[i * n + j];
        Real v = f.v[i * n + j];
        Real norm = std::sqrt(u * u + v * v);
        velocity_sum += norm;
      }
    }
    tr->info() << "VelocitySum=" << velocity_sum;
    Real ref_velocity_sum = 32197.2860476546;
    bool do_compare = true;
    if (do_compare) {
      if (!math::isNearlyEqual(velocity_sum, ref_velocity_sum))
        ARCANE_FATAL("Bad reference ref={0} v={1}", ref_velocity_sum, velocity_sum);
    }
  }
  else
    tr->info() << "Skip velocity test because number of time step or number of cells are not equal to the reference";
}

int
main(int argc,char* argv[])
{
  // Le nom du fichier du jeu de données est le dernier argument de la ligne de commande.
  if (argc<2){
    std::cout << "Usage: SimpleEulerianCFD casefile.arc\n";
    return 1;
  }
  ArcaneLauncher::init(CommandLineArguments(&argc,&argv));
  String case_file_name = argv[argc-1];
  // Déclare la fonction qui sera exécutée par l'appel à run()
  auto f = [=](DirectSubDomainExecutionContext& ctx) -> int
  {
    executeCode(ctx.subDomain());
    return 0;
  };
  // Exécute le fonctor 'f'.
  return ArcaneLauncher::run(f);
}
