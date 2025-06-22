// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "AcceleratorUtils.h"

#include <arcane_version.h>
#include <arcane/utils/IMemoryRessourceMng.h>
#if ARCANE_VERSION>31500
#include <arcane/utils/MemoryUtils.h>
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IMemoryAllocator* AcceleratorUtils::
getMemoryAllocator(eMemoryRessource x)
{
#if ARCANE_VERSION>31500
  return MemoryUtils::getAllocator(x);
#else
  return platform::getDataMemoryRessourceMng()->getAllocator(x);
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

IMemoryAllocator* AcceleratorUtils::
getAcceleratorHostMemoryAllocator()
{
#if ARCANE_VERSION>31500
  return MemoryUtils::getDefaultDataAllocator();
#else
  return platform::getAcceleratorHostMemoryAllocator();
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

