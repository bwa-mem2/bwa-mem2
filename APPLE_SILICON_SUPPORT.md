# Apple Silicon Support for BWA-MEM2

## Modifications apportées

Cette branche ajoute le support complet pour les processeurs Apple Silicon (ARM64) sur macOS.

### Changements principaux

1. **Ajout de sse2neon** : Bibliothèque de traduction SSE vers NEON
   - Ajoutée comme submodule git dans `ext/sse2neon`
   - Permet d'émuler les instructions SIMD x86 sur ARM

2. **Nouveau fichier `src/simd_compat.h`** : Wrapper de compatibilité SIMD
   - Détecte automatiquement l'architecture (x86 vs ARM)
   - Inclut les headers appropriés selon la plateforme
   - Définit des macros de compatibilité (`__rdtsc`, `_mm_prefetch`, etc.)

3. **Modifications du Makefile**
   - Détection automatique de l'architecture ARM64
   - Compilation avec `-march=native` pour ARM
   - Ajout du binaire `bwa-mem2.arm64`

4. **Modifications des fichiers source**
   - `src/fastmap.cpp` : Gestion ARM pour `cpuid` et `HTStatus()`
   - `src/bwamem.cpp` : Code conditionnel pour AVX-512
   - `src/runsimd.cpp` : Lancement du binaire `.arm64` sur ARM
   - Remplacement des includes x86 par `simd_compat.h`

5. **Corrections de compatibilité**
   - `ext/safestringlib` : Fix pour `memset_s` et `abort()` sur macOS
   - Ajout de `.DS_Store` au `.gitignore`

## Compilation

Sur macOS avec Apple Silicon :
```bash
make clean
make
```

Le système détectera automatiquement ARM64 et compilera le binaire approprié.

## Test

```bash
./bwa-mem2
```

Devrait afficher :
```
Looking to launch executable "/path/to/bwa-mem2.arm64", simd = .arm64
Launching executable "/path/to/bwa-mem2.arm64"
Usage: bwa-mem2 <command> <arguments>
```

## Compatibilité

- ✅ macOS avec Apple Silicon (M1/M2/M3/M4)
- ✅ Rétrocompatible avec x86_64
- ✅ Fonctionne sur Intel Mac et Linux x86_64

## Performance

Les instructions NEON d'ARM sont équivalentes en performance aux instructions SSE, donc les performances devraient être similaires ou meilleures grâce à l'architecture Apple Silicon.
