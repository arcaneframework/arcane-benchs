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

constexpr Int32 U_FIELD = 0;
constexpr Int32 V_FIELD = 1;
constexpr Int32 S_FIELD = 2;

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
  static void _destroy(Real** ptr)
  {
    delete *ptr;
    *ptr = nullptr;
  }
  static void _allocateAndFillZero(NumArray<Real, MDDim1>* ptr, Int32 size)
  {
    ptr->resize(size);
    SmallSpan<Real> p = ptr->to1DSmallSpan();
    for (Int32 i = 0; i < size; ++i)
      p[i] = {};
  }
  static void _destroy([[maybe_unused]] NumArray<Real, MDDim1>* ptr)
  {
  }

  static void _doCopy(Real* new_v, const Real* v, Int32 size)
  {
    for (Int32 i = 0; i < size; ++i)
      new_v[i] = v[i];
  }
  static void _doCopy(NumArray<Real, MDDim1>& new_v, const NumArray<Real, MDDim1>& v, [[maybe_unused]] Int32 size)
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

  Real m_dt = 1.0 / 120.0;
  Int32 m_num_iteration = 100;
  Int32 frameNr = 0;
  Real m_over_relaxation = 1.9;
  Real m_obstacle_x = 0.0;
  Real m_obstacle_y = 0.0;
  Real m_obstacle_radius = 0.15;
  Int32 sceneNr = 0;
  std::unique_ptr<Fluid> m_fluid;

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
    m_density = density_;
    m_nb_cell_x = numX_ + 2;
    m_nb_cell_y = numY_ + 2;
    m_height = h_;
    scene = scene_;

    m_nb_cell = m_nb_cell_x * m_nb_cell_y;
    Utils::_allocateAndFillZero(&m_velocity_u, m_nb_cell);
    Utils::_allocateAndFillZero(&m_velocity_v, m_nb_cell);
    Utils::_allocateAndFillZero(&m_new_velocity_u, m_nb_cell);
    Utils::_allocateAndFillZero(&m_new_velocity_v, m_nb_cell);
    Utils::_allocateAndFillZero(&m_p, m_nb_cell);
    Utils::_allocateAndFillZero(&m_is_fluid, m_nb_cell);
    Utils::_allocateAndFillZero(&m_m, m_nb_cell);
    Utils::_allocateAndFillZero(&m_new_m, m_nb_cell);
  }

  ~Fluid()
  {
    Utils::_destroy(&m_velocity_u);
    Utils::_destroy(&m_velocity_v);
    Utils::_destroy(&m_new_velocity_u);
    Utils::_destroy(&m_new_velocity_v);
    Utils::_destroy(&m_p);
    Utils::_destroy(&m_is_fluid);
    Utils::_destroy(&m_m);
    Utils::_destroy(&m_new_m);

  }
  void integrate(Real dt, Real gravity);
  void solveIncompressibility(Int32 numIters, Real dt);
  void extrapolate();
  Real sampleField(Real x, Real y, RealArrayType& f, Int32 field);

  Real avgU(Int32 i, Int32 j)
  {
    Int32 n = m_nb_cell_y;
    Real lu = (m_velocity_u[i * n + j - 1] + m_velocity_u[i * n + j] + m_velocity_u[(i + 1) * n + j - 1] + m_velocity_u[(i + 1) * n + j]) * 0.25;
    return lu;
  }

  Real avgV(Int32 i, Int32 j)
  {
    Int32 n = m_nb_cell_y;
    Real lv = (m_velocity_v[(i - 1) * n + j] + m_velocity_v[i * n + j] + m_velocity_v[(i - 1) * n + j + 1] + m_velocity_v[i * n + j + 1]) * 0.25;
    return lv;
  }

  void advectVelocity(Real dt);
  void advectSmoke(Real dt);
  void simulate(Real dt, Int32 numIters);

 public:

  Int32 m_nb_cell_x = 0;
  Int32 m_nb_cell_y = 0;
  Int32 m_nb_cell = 0;
  Real m_density = 0.0;
  Real m_height = 0.0;
  RealArrayType m_velocity_u;
  RealArrayType m_velocity_v;
  RealArrayType m_new_velocity_u;
  RealArrayType m_new_velocity_v;
  RealArrayType m_p;
  RealArrayType m_is_fluid;
  RealArrayType m_m;
  RealArrayType m_new_m;
  Scene* scene = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
solveIncompressibility(Int32 numIters, Real dt)
{
  const Int32 n = m_nb_cell_y;
  Real cp = m_density * m_height / dt;

  for (Int32 iter = 0; iter < numIters; iter++) {

    for (Int32 i = 1; i < m_nb_cell_x - 1; i++) {
      for (Int32 j = 1; j < m_nb_cell_y - 1; j++) {

        if (m_is_fluid[i * n + j] == 0.0)
          continue;

        Real sx0 = m_is_fluid[(i - 1) * n + j];
        Real sx1 = m_is_fluid[(i + 1) * n + j];
        Real sy0 = m_is_fluid[i * n + j - 1];
        Real sy1 = m_is_fluid[i * n + j + 1];
        Real ls = sx0 + sx1 + sy0 + sy1;

        if (ls == 0.0)
          continue;

        Real div = m_velocity_u[(i + 1) * n + j] - m_velocity_u[i * n + j] + m_velocity_v[i * n + j + 1] - m_velocity_v[i * n + j];
        Real lp = -div / ls;
        lp *= scene->m_over_relaxation;
        m_p[i * n + j] += cp * lp;
        m_velocity_u[i * n + j] -= sx0 * lp;
        m_velocity_u[(i + 1) * n + j] += sx1 * lp;
        m_velocity_v[i * n + j] -= sy0 * lp;
        m_velocity_v[i * n + j + 1] += sy1 * lp;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
extrapolate()
{
  const Int32 n = m_nb_cell_y;
  for (Int32 i = 0; i < m_nb_cell_x; i++) {
    m_velocity_u[i * n + 0] = m_velocity_u[i * n + 1];
    m_velocity_u[i * n + m_nb_cell_y - 1] = m_velocity_u[i * n + m_nb_cell_y - 2];
  }
  for (Int32 j = 0; j < m_nb_cell_y; j++) {
    m_velocity_v[0 * n + j] = m_velocity_v[1 * n + j];
    m_velocity_v[(m_nb_cell_x - 1) * n + j] = m_velocity_v[(m_nb_cell_x - 2) * n + j];
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Real Fluid::
sampleField(Real x, Real y, RealArrayType& f, Int32 field)
{
  Int32 n = m_nb_cell_y;
  Real h1 = 1.0 / m_height;
  Real h2 = 0.5 * m_height;

  x = std::max(std::min(x, m_nb_cell_x * m_height), m_height);
  y = std::max(std::min(y, m_nb_cell_y * m_height), m_height);

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

  Int32 x0 = std::min(Utils::_doInt32Floor((x - dx) * h1), m_nb_cell_x - 1);
  Real tx = ((x - dx) - x0 * m_height) * h1;
  Int32 x1 = std::min(x0 + 1, m_nb_cell_x - 1);

  Int32 y0 = std::min(Utils::_doInt32Floor((y - dy) * h1), m_nb_cell_y - 1);
  Real ty = ((y - dy) - y0 * m_height) * h1;
  Int32 y1 = std::min(y0 + 1, m_nb_cell_y - 1);

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
advectVelocity(Real dt)
{
  Utils::_doCopy(m_new_velocity_u, m_velocity_u, m_nb_cell);
  Utils::_doCopy(m_new_velocity_v, m_velocity_v, m_nb_cell);

  Int32 n = m_nb_cell_y;
  Real h2 = 0.5 * m_height;

  for (Int32 i = 1; i < m_nb_cell_x; i++) {
    for (Int32 j = 1; j < m_nb_cell_y; j++) {

      // u component
      if (m_is_fluid[i * n + j] != 0.0 && m_is_fluid[(i - 1) * n + j] != 0.0 && j < m_nb_cell_y - 1) {
        Real x = i * m_height;
        Real y = j * m_height + h2;
        Real lu = m_velocity_u[i * n + j];
        Real lv = avgV(i, j);
        x = x - dt * lu;
        y = y - dt * lv;
        Real new_u = sampleField(x, y, m_velocity_u, U_FIELD);
        m_new_velocity_u[i * n + j] = new_u;
      }
      // v component
      if (m_is_fluid[i * n + j] != 0.0 && m_is_fluid[i * n + j - 1] != 0.0 && i < m_nb_cell_x - 1) {
        Real x = i * m_height + h2;
        Real y = j * m_height;
        Real lu = avgU(i, j);
        Real lv = m_velocity_v[i * n + j];
        x = x - dt * lu;
        y = y - dt * lv;
        Real new_v = sampleField(x, y, m_velocity_v, V_FIELD);
        m_new_velocity_v[i * n + j] = new_v;
      }
    }
  }

  Utils::_doCopy(m_velocity_u, m_new_velocity_u, m_nb_cell);
  Utils::_doCopy(m_velocity_v, m_new_velocity_v, m_nb_cell);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
advectSmoke(Real dt)
{
  Utils::_doCopy(m_new_m, m_m, m_nb_cell);

  Int32 n = m_nb_cell_y;
  Real h2 = 0.5 * m_height;

  for (Int32 i = 1; i < m_nb_cell_x - 1; i++) {
    for (Int32 j = 1; j < m_nb_cell_y - 1; j++) {

      if (m_is_fluid[i * n + j] != 0.0) {
        Real lu = (m_velocity_u[i * n + j] + m_velocity_u[(i + 1) * n + j]) * 0.5;
        Real lv = (m_velocity_v[i * n + j] + m_velocity_v[i * n + j + 1]) * 0.5;
        Real x = i * m_height + h2 - dt * lu;
        Real y = j * m_height + h2 - dt * lv;

        m_new_m[i * n + j] = sampleField(x, y, m_m, S_FIELD);
      }
    }
  }
  Utils::_doCopy(m_m, m_new_m, m_nb_cell);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Fluid::
simulate(Real dt, Int32 numIters)
{
  Utils::_doFillZero(m_p, m_nb_cell);
  solveIncompressibility(numIters, dt);

  extrapolate();
  advectVelocity(dt);
  advectSmoke(dt);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Scene::
setupScene(Int32 resolution)
{
  m_obstacle_radius = 0.15;
  m_over_relaxation = 1.9;

  m_dt = 1.0 / 60.0;
  m_num_iteration = 40;

  Real domainHeight = 1.0;
  Real domainWidth = domainHeight / simHeight * simWidth;
  Real h = domainHeight / resolution;

  Int32 numX = Utils::_doInt32Floor(domainWidth / h);
  Int32 numY = Utils::_doInt32Floor(domainHeight / h);
  std::cout << "domain_height=" << domainHeight << " width=" << domainWidth << " computed_h=" << h << "\n";

  Real density = 1000.0;

  std::cout << "NumX=" << numX << " NumY=" << numY << " h=" << h << "\n";

  m_fluid = std::make_unique<Fluid>(density, numX, numY, h, this);
  Fluid& f = *m_fluid;

  const Int32 n = f.m_nb_cell_y;

  {
    // vortex shedding
    Real inVel = 2.0;
    for (Int32 i = 0; i < f.m_nb_cell_x; i++) {
      for (Int32 j = 0; j < f.m_nb_cell_y; j++) {
        Real s = 1.0; // fluid
        if (i == 0 || j == 0 || j == f.m_nb_cell_y - 1)
          s = 0.0; // solid
        f.m_is_fluid[i * n + j] = s;

        if (i == 1) {
          f.m_velocity_u[i * n + j] = inVel;
        }
      }
    }

    Real pipeH = 0.1 * f.m_nb_cell_y;
    Int32 minJ = Utils::_doInt32Floor(0.5 * f.m_nb_cell_y - 0.5 * pipeH);
    Int32 maxJ = Utils::_doInt32Floor(0.5 * f.m_nb_cell_y + 0.5 * pipeH);

    for (Int32 j = minJ; j < maxJ; j++)
      f.m_m[j] = 0.0;

    setObstacle(0.4, 0.5, true);

    {
      // Vortex shedding
      m_dt = 1.0 / 120.0;
      m_num_iteration = 100;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Scene::
setObstacle(Real x, Real y, bool reset)
{
  Real vx = 0.0;
  Real vy = 0.0;

  if (!reset) {
    vx = (x - m_obstacle_x) / m_dt;
    vy = (y - m_obstacle_y) / m_dt;
  }

  m_obstacle_x = x;
  m_obstacle_y = y;
  Real r = m_obstacle_radius;
  Fluid& f = *m_fluid;
  Int32 n = f.m_nb_cell_y;

  for (Int32 i = 1; i < f.m_nb_cell_x - 2; i++) {
    for (Int32 j = 1; j < f.m_nb_cell_y - 2; j++) {

      f.m_is_fluid[i * n + j] = 1.0;

      Real dx = (i + 0.5) * f.m_height - x;
      Real dy = (j + 0.5) * f.m_height - y;

      if (dx * dx + dy * dy < r * r) {
        f.m_is_fluid[i * n + j] = 0.0;
        f.m_m[i * n + j] = 1.0;
        f.m_velocity_u[i * n + j] = vx;
        f.m_velocity_u[(i + 1) * n + j] = vx;
        f.m_velocity_v[i * n + j] = vy;
        f.m_velocity_v[i * n + j + 1] = vy;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void executeCode(ISubDomain* sd)
{
  // This code does not use a mesh

  // Get the instance to manage listing.
  ITraceMng* tr = sd->traceMng();

  // Number of cells by dimension
  const Int32 resolution = 100;
  // Number of time step
  const Int32 nb_time_step = 100;

  Scene scene;
  scene.setupScene(resolution);
  Real begin_x = Platform::getRealTime();
  for (Int32 i = 0; i < 100; ++i) {
    Fluid& f = *scene.m_fluid;
    if (i < 10 || ((i % 5) == 0))
      tr->info() << "DoSimulate nx=" << f.m_nb_cell_x << " ny=" << f.m_nb_cell_y << " dt=" << scene.m_dt
                 << " frame=" << scene.frameNr;
    f.simulate(scene.m_dt, scene.m_num_iteration);
    ++scene.frameNr;
  }
  Real end_x = Platform::getRealTime();
  tr->info() << "ComputeTime=" << (end_x - begin_x);

  // Check velocity
  // Only valid if resolution == 100 and nb_time_step == 100
  if (nb_time_step == 100 && resolution == 100) {
    Real velocity_sum = 0.0;
    Fluid& f = *scene.m_fluid;
    Int32 n = f.m_nb_cell_y;
    for (Int32 i = 1; i < f.m_nb_cell_x - 2; i++) {
      for (Int32 j = 1; j < f.m_nb_cell_y - 2; j++) {
        Real u = f.m_velocity_u[i * n + j];
        Real v = f.m_velocity_v[i * n + j];
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
