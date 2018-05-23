#ifndef CPPQT_H

#define CPPQT_H

#include <QtCore/QMap>
#include <QtCore/QSet>
#include <QtWidgets/QGraphicsPolygonItem>
#include <QtWidgets/QGraphicsScene>

QString GenerateUuid();

class TriangleInfo;
class PointGraphicsItemBase;
class TriangleItem;

class PointInfo
{
public:
	PointInfo(const QString &uuid, double x = 0, double y = 0, const QString &label = "", double timedelay = qQNaN())
		: m_uuid(uuid), m_x(x), m_y(y), m_label(label), m_timedelay(timedelay) { }
	~PointInfo();

	bool isTimedelayValid() const { return !qIsNaN(m_timedelay); }

	QString m_uuid;
	double m_x, m_y;
	QString m_label;
	double m_timedelay;

	QSet<TriangleInfo *> m_triangles;
};

class TriangleInfo
{
public:
	TriangleInfo(const QString &triangUuid, PointInfo *p0, PointInfo *p1, PointInfo *p2);
	~TriangleInfo();

	QString m_triangleUuid;
	PointInfo *m_p[3];
};

class Layer
{
public:
	Layer(const QString &name = "Untitled");
	virtual ~Layer();

	void setName(const QString &n) { m_name = n; }
	QString getName() const { return m_name; }
	QString getUuid() const { return m_uuid; }

	static QString generateKey() { return GenerateUuid(); }
	QString addPoint(double x, double y, QString label = QString(), double timedelay = qQNaN());
	QString addPoint(QVariant xy, QVariant label = QVariant(), QVariant timedelay = QVariant());
	PointInfo *setPoint(QString uuid, double x, double y, QString label = QString(), double timedelay = qQNaN());
	void setPoint(QString uuid, QVariant xy, QVariant label = QVariant(), QVariant timedelay = QVariant());
	QVariant clearPoint(QString uuid);
	void clearAllPoints();

	PointInfo *getPointInfo(QString uuid);
    const QMap<QString, PointInfo *> &getPointInfos() const { return m_points; }

	QVariant getPoint(QString uuid) const;
	QVariant getPoints(bool transformed = false);

	QTransform getImageTransform() const { return m_transform; }
	void setImageTransform(const QTransform &t) { m_transform = t; }

	QVariant getTrianglesContainingPoint(QString uuid);
	QString addTriangle(const QString &p0Key, const QString &p1Key, const QString &p2Key, const QString &triangUuid = QString());
	TriangleInfo *addTriangleInfo(const QString &p0Key, const QString &p1Key, const QString &p2Key, const QString &triangUuid = QString());
	const TriangleInfo *getTriangleInfo(const QString &triangUuid) const;
	const QMap<QString, TriangleInfo *> &getTriangleInfos() const { return m_triangles; } 
	QVariant getTriangle(const QString &triangUuid) const;
	QVariant getTriangles() const;
	void clearTriangle(const QString &triangUuid);
private:
	QString m_name;
	QTransform m_transform;
	QString m_uuid;
	QMap<QString, PointInfo *> m_points;
	QMap<QString, TriangleInfo *> m_triangles;
};

class EmptyGraphicsItem : public QGraphicsItem
{
public:
	EmptyGraphicsItem(bool childrenBoundingRect, QGraphicsItem *pParent = 0);
	~EmptyGraphicsItem();

	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
private:
	bool m_childrenBoundingRect;
};

class TriangleItem : public QGraphicsItem
{
public:
	TriangleItem(QGraphicsItem *pParent = 0);
	TriangleItem(const QString &uuid, Layer *pLayer, const QMap<QString, QPointF> &pointCoords, QGraphicsItem *pParent = 0);
	~TriangleItem();

	QString getUuid() const { return m_triangleUuid; }
	Layer *getLayer() { return m_pLayer; }
	void updatePointPosition(const QString &ptUuid, QPointF &pos);

	bool contains(const QPointF &p) const;
	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);

	static TriangleItem *getTriangleItem(QGraphicsItem *pItem);

	void clearMarker() { m_marker = false; }
	void mark() { m_marker = true; }
	bool isMarked() const { return m_marker; }
private:
	void updatePolygon();

	static QBrush s_normalBrush, s_selectedBrush;
	static QPen s_normalPen, s_selectedPen;
	static bool s_brushInit;

	QPolygonF m_polygon;
	QString m_triangleUuid;
	Layer *m_pLayer;
	QMap<QString, QPointF> m_ptsCoords;

	bool m_marker;
};

class SyncTextEditItem : public QGraphicsTextItem
{
public:
	SyncTextEditItem(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent);
	~SyncTextEditItem();

	void recenter();
private:
	void keyReleaseEvent(QKeyEvent *event);
	void inputMethodEvent(QInputMethodEvent *event);
	void focusOutEvent(QFocusEvent *event);

	void syncValueAndCenter();

	Layer *m_pLayer;
	QString m_pointUuid;
};

class LayerObjectGraphicsItem : public EmptyGraphicsItem
{
public:
	LayerObjectGraphicsItem();
	LayerObjectGraphicsItem(Layer *pLayer, const QString &objectUuid, QGraphicsItem *pParent);
	~LayerObjectGraphicsItem();

	const QString &getUuid() const { return m_uuid; }
	Layer *getLayer() { return m_pLayer; }
private:
	Layer *m_pLayer;
	QString m_uuid;
};

class PointGraphicsItemBase : public LayerObjectGraphicsItem
{
public:
	PointGraphicsItemBase();
	PointGraphicsItemBase(Layer *pLayer, const QString &uuid, const QString &label, double timedelay, QGraphicsItem *pParent);
	~PointGraphicsItemBase();

	bool isMovable() const;
	bool isNormalPoint() const;
	void setMovable(bool v = true);
	void fetchSettings();
	void fetchSettings(const PointInfo &ptInf);
	QVariant syncPosition();

	virtual void onSelected(bool v) { }
	QVariant itemChange(GraphicsItemChange change, const QVariant &value);

	void toggleFocus();
	
	static PointGraphicsItemBase *getPointGraphicsItem(QGraphicsItem *pItem);

	void clearMarker() { m_marker = false; }
	void mark() { m_marker = true; }
	bool isMarked() const { return m_marker; }
private:
	void setTDText(double timedelay);

	static QBrush s_fontBrush;

	SyncTextEditItem *m_pTxt;
	QGraphicsSimpleTextItem *m_pTdTxt;
	bool m_marker;
};

class SinglePointGraphicsItem : public PointGraphicsItemBase
{
public:
	SinglePointGraphicsItem(Layer *pLayer, const QString &uuid, const QString &label, double timedelay, QGraphicsItem *pParent);
	~SinglePointGraphicsItem();

	void onSelected(bool v);

	QRectF boundingRect() const;
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);
private:
	static QBrush s_normalBrush, s_selectedBrush;
	static QPen s_normalPen, s_selectedPen;

	QGraphicsEllipseItem *m_pCircle;
	QGraphicsLineItem *m_pLine1, *m_pLine2;
};

class MatchPointGraphicsItemBase : public PointGraphicsItemBase
{
public:
	MatchPointGraphicsItemBase(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent);
	~MatchPointGraphicsItemBase();
protected:
	static QPen s_normalLinePen, s_selectedLinePen, s_rectPen;
	static QBrush s_rectBrush;
private:
	QGraphicsRectItem *m_pRect;
};

class MatchPointGraphicsItemCross : public MatchPointGraphicsItemBase
{
public:
	MatchPointGraphicsItemCross(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent);
	~MatchPointGraphicsItemCross();

	void onSelected(bool v);
private:
	QGraphicsLineItem *m_pLine1, *m_pLine2;
};

class MatchPointGraphicsItemCircle : public MatchPointGraphicsItemBase
{
public:
	MatchPointGraphicsItemCircle(Layer *pLayer, const QString &uuid, const QString &label, QGraphicsItem *pParent);
	~MatchPointGraphicsItemCircle();

	void onSelected(bool v);
private:
	QGraphicsEllipseItem *m_pCircle;
};

class SceneBase : public QGraphicsScene
{
public:
	SceneBase() { }
	~SceneBase() { }

	void setPointTransform(const QTransform &t);
	virtual void onPointLabelChanged(const QString &layerUuid, const QString &ptUuid, const QString &oldLabel, const QString &newLabel) { }
	
	QTransform getPointTransform() const { return m_pointTransform; }
private:
	QTransform m_pointTransform;
};

class LayerGraphicsItemBase : public EmptyGraphicsItem
{
public:
	enum PointType { Normal, Circle, Cross };

	LayerGraphicsItemBase(Layer *pLayer, PointType pType, bool childrenBoundingRect, QGraphicsItem *pParent = 0);
	~LayerGraphicsItemBase();

	Layer *getLayer() { return m_pLayer; }

	void onPointPositionChanged(const QString &pointUuid, QPointF pos);
	void setPointsVisible(bool v);
	void updatePoints();

    void updatePointList(const QList<QVariant> &points, const QList<QVariant> &triangles);
    void setLabel(const QString &uuid, const QString &label);
    void movePoint(const QString &uuid, double x, double y);
    QString addPoint(double x, double y, const QString &label = QString(), double timedelay = qQNaN(), const QString &uuid = QString());
    void setPoint(const QString &uuid, double x, double y, const QString &label, double timedelay);
    QString addTriangle(const QString &p0, const QString &p1, const QString &p2, const QString &triangleUuid = QString());

    void clearTriangle(const QString &triangleUuid);
    QVariant clearPoint(const QString &uuid);
protected:
	PointGraphicsItemBase *createPoint(const PointInfo &pt);
	TriangleItem *createTriangle(const TriangleInfo *triangInf);

	Layer *m_pLayer;
	PointType m_pointType;
    bool m_pointsVisible;

	QMap<QString, PointGraphicsItemBase *> m_pointItems;
	QMap<QString, TriangleItem *> m_triangleItems;
};

#endif // CPPQT_H

