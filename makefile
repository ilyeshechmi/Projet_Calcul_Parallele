# ===========================
#     Makefile PROPRE
# ===========================

FC      = gfortran
FFLAGS  = -O2 -Wall -std=f2008

# Répertoires
SRC_DIR   = src
BUILD_DIR = build

TARGET  = run

# ===========================
#  Sources (ordre logique)
# ===========================
SRC_FILES = \
	precision.f90 \
	initialiser.f90 \
	donnees.f90 \
	euler.f90 \
	fonctions.f90 \
	exact_sod.f90 \
	schema_rusanov.f90 \
	main.f90

SRCS = $(addprefix $(SRC_DIR)/,$(SRC_FILES))
OBJS = $(SRCS:$(SRC_DIR)/%.f90=$(BUILD_DIR)/%.o)

# ===========================
#       Cible principale
# ===========================
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# ===========================
#   Compilation des .o
# ===========================
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J $(BUILD_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# ===========================
#   Dépendances modules
# ===========================
$(BUILD_DIR)/initialiser.o: $(BUILD_DIR)/precision.o
$(BUILD_DIR)/donnees.o:     $(BUILD_DIR)/precision.o

$(BUILD_DIR)/euler.o:       $(BUILD_DIR)/precision.o
$(BUILD_DIR)/fonctions.o:   $(BUILD_DIR)/precision.o $(BUILD_DIR)/euler.o
$(BUILD_DIR)/exact_sod.o:   $(BUILD_DIR)/precision.o $(BUILD_DIR)/euler.o

$(BUILD_DIR)/schema_rusanov.o: \
	$(BUILD_DIR)/precision.o \
	$(BUILD_DIR)/fonctions.o \
	$(BUILD_DIR)/euler.o

$(BUILD_DIR)/main.o: \
	$(BUILD_DIR)/precision.o \
	$(BUILD_DIR)/initialiser.o \
	$(BUILD_DIR)/donnees.o \
	$(BUILD_DIR)/fonctions.o \
	$(BUILD_DIR)/schema_rusanov.o \
	$(BUILD_DIR)/euler.o \
	$(BUILD_DIR)/exact_sod.o

# ===========================
#       Nettoyage
# ===========================
clean:
	rm -rf $(BUILD_DIR) *.dat $(TARGET) resultats/*.dat

.PHONY: all clean
