%Module cppqt

%Import QtCore/QtCoremod.sip
%Import QtWidgets/QtWidgetsmod.sip

class Layer
{
%TypeHeaderCode
#include "cppqt.h"
%End

public:
	Layer(const QString &name = "Untitled");
	virtual ~Layer();

	void setName(const QString &n);
	QString getName() const;
	QString getUuid() const;

	static QString generateKey();
	QString addPoint(QVariant xy, QVariant label = QVariant(), QVariant timedelay = QVariant());
	void setPoint(QString uuid, QVariant xy, QVariant label = QVariant(), QVariant timedelay = QVariant());
	QVariant clearPoint(QString uuid);

	QVariant getPoint(QString uuid) const;
	QVariant getPoints(bool transformed = false);

	QTransform getImageTransform() const;
	void setImageTransform(const QTransform &t);

	QVariant getTrianglesContainingPoint(QString uuid);
	QString addTriangle(const QString &p0Key, const QString &p1Key, const QString &p2Key, const QString &triangUuid = QString());
	QVariant getTriangle(const QString &triangUuid) const;
	QVariant getTriangles() const;
	void clearTriangle(const QString &triangUuid);
};

class EmptyGraphicsItem : public QGraphicsItem
{

%TypeHeaderCode
#include "cppqt.h"
%End

public:
	EmptyGraphicsItem(bool childrenBoundingRect, QGraphicsItem *pParent /TransferThis/ = 0);
	~EmptyGraphicsItem();

	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
};

class TriangleItem : public QGraphicsItem
{
%TypeHeaderCode
#include "cppqt.h"
%End

public:
	TriangleItem(QGraphicsItem *pParent /TransferThis/ = 0);
	~TriangleItem();

	QString getUuid() const;
	Layer *getLayer();

	static TriangleItem *getTriangleItem(QGraphicsItem *pItem);
};

class LayerObjectGraphicsItem : public EmptyGraphicsItem
{
%TypeHeaderCode
#include "cppqt.h"
%End
public:
	LayerObjectGraphicsItem();
	~LayerObjectGraphicsItem();

	const QString &getUuid() const;
	Layer *getLayer();
};

class PointGraphicsItemBase : public LayerObjectGraphicsItem
{
%TypeHeaderCode
#include "cppqt.h"
%End
public:
	PointGraphicsItemBase();
	~PointGraphicsItemBase();

	bool isMovable();
	bool isNormalPoint();
	void setMovable(bool v = true);
	void fetchSettings();
	QVariant syncPosition();
	void toggleFocus();

	static PointGraphicsItemBase *getPointGraphicsItem(QGraphicsItem *pItem);
};

class SceneBase : public QGraphicsScene
{
%TypeHeaderCode
#include "cppqt.h"
%End
public:
	SceneBase();
	~SceneBase();

	void setPointTransform(const QTransform &t);
	QTransform getPointTransform() const;

	virtual void onPointLabelChanged(const QString &layerUuid, const QString &ptUuid, const QString &oldLabel, const QString &newLabel);
};

class LayerGraphicsItemBase : public EmptyGraphicsItem
{
%TypeHeaderCode
#include "cppqt.h"
%End
public:
	enum PointType { Normal, Circle, Cross };

	LayerGraphicsItemBase(Layer *pLayer, PointType pType, bool childrenBoundingRect, QGraphicsItem *pParent /TransferThis/ = 0);
	~LayerGraphicsItemBase();

	Layer *getLayer();
	void updatePoints();
	void setPointsVisible(bool v);

	void updatePointList(const QList<QVariant> &points, const QList<QVariant> &triangles);
	void setLabel(const QString &uuid, const QString &label);
	void movePoint(const QString &uuid, double x, double y);
	QString addPoint(double x, double y, const QString &label = QString(), double timedelay = qQNaN(), const QString &uuid = QString());
	void setPoint(const QString &uuid, double x, double y, const QString &label, double timedelay);
	QString addTriangle(const QString &p0, const QString &p1, const QString &p2, const QString &triangleUuid = QString());

	void clearTriangle(const QString &triangleUuid);
	QVariant clearPoint(const QString &uuid);
};

%ModuleCode
#include "cppqt.h"
#include "cppqt.cpp"
%End
