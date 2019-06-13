#ifndef NEW_NET_GLOBAL_H
#define NEW_NET_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(NEW_NET_LIBRARY)
#  define NEW_NETSHARED_EXPORT Q_DECL_EXPORT
#else
#  define NEW_NETSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // NEW_NET_GLOBAL_H
